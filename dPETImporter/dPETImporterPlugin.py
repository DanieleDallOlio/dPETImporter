import re, logging
import vtk, qt, ctk, slicer
import vtk.util.numpy_support
import DICOMLib
from DICOMLib import DICOMPlugin, DICOMLoadable
from slicer.util import settingsValue, toBool
import json
from datetime import datetime, timedelta

#
# This is the plugin to efficiently handle translation of DICOM objects
# that represent 2D or 3D dynamic PET acquisitions. This plugin converts such
# DICOM files into MRML nodes as multivolume objects. We follow the
# MultiVolumeImporterPlugin as reference.
#

def _parse_dicom_datetime(date_str, time_str):
  if not date_str or not time_str:
    return None
  for fmt in ("%Y%m%d%H%M%S.%f", "%Y%m%d%H%M%S"):
    try:
      return datetime.strptime(date_str + time_str, fmt)
    except Exception:
      pass
  return None


def _parse_dicom_datetime_value(dt_str):
  if not dt_str:
    return None
  text = str(dt_str).strip()
  # DICOM DT may contain fractional seconds and an optional UTC offset.
  # For the injection-to-acquisition interval we need the local wall-clock
  # value; acquisition DT in many PET objects does not carry an offset, so
  # strip the optional suffix consistently rather than mixing aware/naive DTs.
  m = re.match(r"^(\d{14})(?:\.(\d+))?(?:[+-]\d{4})?$", text)
  if not m:
    return None
  base = m.group(1)
  frac = m.group(2)
  try:
    dt = datetime.strptime(base, "%Y%m%d%H%M%S")
    if frac:
      microseconds = int((frac[:6]).ljust(6, '0'))
      dt = dt.replace(microsecond=microseconds)
    return dt
  except Exception:
    return None


def compute_suvbw_factor(mvNode):
  """
  Compute the Bq/mL -> SUVbw factor for the DICOM decay convention.

  Supported conventions:
    START: dose is decayed from administration to the first acquisition start.
    ADMIN: administered dose is used directly; administration datetime is NOT
           used for any back/forward extrapolation.

  NONE/unknown are intentionally unsupported for quantitative Multi-timepoint
  use and return None.
  """
  try:
    half_life = float(mvNode.GetAttribute('RadionuclideHalfLife') or 0.0)
    weight = float(mvNode.GetAttribute('PatientWeight') or 0.0)
    totalDose = float(mvNode.GetAttribute('RadionuclideTotalDose') or 0.0)
    decayCorrection = (mvNode.GetAttribute('DecayCorrection') or '').upper()
    if weight <= 0.0 or totalDose <= 0.0:
      return None
    dose_kbq = totalDose * 0.001
    if dose_kbq <= 0.0:
      return None

    if decayCorrection == 'ADMIN':
      # ADMIN-referenced image activity is normalized by the administered
      # (non-decayed) dose.  No administration datetime is used here.
      return float(weight / dose_kbq)

    if decayCorrection != 'START':
      return None

    rstart = mvNode.GetAttribute('RadionuclideStartDateTime')
    first_frame_dt = mvNode.GetAttribute('dPET.FirstFrameAcquisitionDateTime')
    if not (rstart and first_frame_dt and half_life > 0.0):
      return None

    start_dt = _parse_dicom_datetime_value(rstart)
    ref_dt = _parse_dicom_datetime_value(first_frame_dt)
    if start_dt is None or ref_dt is None:
      return None

    decay_seconds = (ref_dt - start_dt).total_seconds()
    decayedDose = dose_kbq * 2 ** (-decay_seconds / half_life)
    if decayedDose <= 0.0:
      return None
    return float(weight / decayedDose)
  except Exception:
    return None


def compute_suvbw_for_start(mvNode):
  """Backward-compatible wrapper retained for external callers."""
  if (mvNode.GetAttribute('DecayCorrection') or '').upper() != 'START':
    return None
  return compute_suvbw_factor(mvNode)



# PET volumes can have a very sparse intensity histogram (notably early frames
# of dynamic PET).  Slicer's 0.1/99.9-percentile automatic window can then
# collapse even though the image contains valid nonzero voxel values.  Keep the
# ordinary automatic behavior unless that exceptional condition is detected.
def _dpetApplyCollapsedAutoWindowFallback(volumeNode):
  """Fix only PET auto W/L failures caused by empty or extremely sparse frames.

  Voxel data are never modified. A sparse frame means the top 0.1 percent
  excluded by Slicer's histogram W/L contain essentially all meaningful
  activity. The frame-specific display setting is managed by the sequence
  browser observer (which restores auto W/L on the next normal frame).
  """
  import math
  import numpy as np
  if volumeNode is None:
    return None
  display = volumeNode.GetDisplayNode()
  image = volumeNode.GetImageData()
  if display is None or image is None or not display.GetAutoWindowLevel():
    return None

  valueType = (volumeNode.GetAttribute('dPET.ValueType') or '').upper()
  if valueType == 'SUVBW':
    displayFloor = 1.0
  elif valueType == 'BQML':
    displayFloor = 1000.0
  else:
    return None

  # CalculateAutoLevels() is a C++ MRML method that is not wrapped in Slicer's
  # Python API. Modified() invokes Slicer's native auto-W/L recalculation
  # for the current volume while AutoWindowLevel is enabled.
  display.Modified()
  window = float(display.GetWindow())
  level = float(display.GetLevel())
  low, high = (float(x) for x in image.GetScalarRange())
  if not all(math.isfinite(v) for v in (low, high, window, level)):
    return None

  reason = None
  activeFraction = None
  autoUpper = level + window / 2.0
  if high <= displayFloor * 1.0e-6:
    # An almost-empty volume must NOT map quantization noise to white.
    reason = 'near-zero image'
  elif window <= displayFloor * 1.0e-8:
    reason = 'numerically collapsed window'
  elif (window < displayFloor and autoUpper < displayFloor):
    # An apparently nonzero auto window can still show only numerical noise
    # when <0.1% of voxels carry measurable PET activity. This is common in
    # early dynamic PET where Slicer's 99.9th percentile misses the signal.
    scalars = image.GetPointData().GetScalars() if image.GetPointData() else None
    if scalars is not None and scalars.GetNumberOfComponents() == 1:
      samples = vtk.util.numpy_support.vtk_to_numpy(scalars)
      stride = max(1, samples.size // 250000)
      samples = samples[::stride]
      samples = samples[np.isfinite(samples)]
      if samples.size:
        signalThreshold = min(0.01 * displayFloor, 1.0e-4 * high)
        activeFraction = float(np.count_nonzero(samples > signalThreshold) / samples.size)
        if activeFraction < 0.001:
          reason = 'meaningful activity excluded from auto histogram'

  if reason is None:
    return None

  upper = max(displayFloor, high)
  if not math.isfinite(upper) or upper <= 0.0:
    return None
  display.AutoWindowLevelOff()
  display.SetWindowLevelMinMax(0.0, upper)
  logging.info('[dPET display] PET W/L fallback: %s, autoW=%.10g, '
               'autoL=%.10g, image=[%.10g, %.10g], signalFraction=%s; '
               'fallback=[0, %.10g] %s', reason, window, level, low, high,
               activeFraction, upper, valueType)
  return (float(display.GetWindow()), float(display.GetLevel()))


_dpetPETWindowGuards = {}


class _dpetPETSequenceWindowGuard:
  """Recheck automatic PET windowing after each Sequence Browser item change.

  The display node belongs to the proxy, not to each stored frame.  Therefore
  a fallback must be undone before evaluating the next item.  A user-entered
  manual window/level always takes precedence.
  """
  def __init__(self, browser, sequence, proxy):
    self.browser = browser
    self.sequence = sequence
    self.proxy = proxy
    self.browserID = browser.GetID()
    self.lastItem = None
    self.fallbackWL = None
    self.queued = False
    self.active = True
    self.browserObserver = browser.AddObserver(
      vtk.vtkCommand.ModifiedEvent, self._onBrowserModified)
    self.browserDeleteObserver = browser.AddObserver(
      vtk.vtkCommand.DeleteEvent, self._onBrowserDeleted)
    self.sceneObserver = slicer.mrmlScene.AddObserver(
      slicer.vtkMRMLScene.StartCloseEvent, self._onSceneClosing)
    self._schedule()

  def _onBrowserModified(self, caller=None, event=None):
    if self.active and self.browser.GetSelectedItemNumber() != self.lastItem:
      self._schedule()

  def _onSceneClosing(self, caller=None, event=None):
    self.close()

  def _onBrowserDeleted(self, caller=None, event=None):
    # Do not touch observer tags on an object being destroyed.
    self.close(browserDeleted=True)

  def _schedule(self):
    if self.queued or not self.active:
      return
    self.queued = True
    # Slicer's sequence logic must finish copying the new frame into the proxy
    # before its automatic W/L is inspected.
    qt.QTimer.singleShot(0, self._evaluate)

  @staticmethod
  def _sameWindowLevel(display, expected):
    return (abs(float(display.GetWindow()) - expected[0]) <=
            1.0e-9 * max(1.0, abs(expected[0]))
            and abs(float(display.GetLevel()) - expected[1]) <=
            1.0e-9 * max(1.0, abs(expected[1])))

  def _evaluate(self):
    self.queued = False
    if not self.active:
      return
    try:
      if (self.browser.GetScene() is None or self.proxy.GetScene() is None
          or self.browser.GetMasterSequenceNode() != self.sequence):
        self.close()
        return
      item = self.browser.GetSelectedItemNumber()
      if item < 0 or item == self.lastItem:
        return
      display = self.proxy.GetDisplayNode()
      if display is None:
        return

      if self.fallbackWL is not None:
        # Restore automatic W/L only if the current manual window is exactly
        # the one installed by this guard.  If the user adjusted it, preserve
        # their manual setting across subsequent frames.
        if (not display.GetAutoWindowLevel()
            and self._sameWindowLevel(display, self.fallbackWL)):
          display.AutoWindowLevelOn()
        self.fallbackWL = None

      if display.GetAutoWindowLevel():
        self.fallbackWL = _dpetApplyCollapsedAutoWindowFallback(self.proxy)
      # Otherwise the user has selected manual W/L.  Never change it.
      self.lastItem = item
    except Exception:
      logging.warning('[dPET display] Conditional PET W/L safeguard failed',
                      exc_info=True)

  def close(self, browserDeleted=False):
    if not self.active:
      return
    self.active = False
    _dpetPETWindowGuards.pop(self.browserID, None)
    if not browserDeleted:
      try:
        self.browser.RemoveObserver(self.browserObserver)
        self.browser.RemoveObserver(self.browserDeleteObserver)
      except Exception:
        pass
    try:
      slicer.mrmlScene.RemoveObserver(self.sceneObserver)
    except Exception:
      pass


def _dpetWatchPETSequenceDisplay(browser, sequence, proxy):
  """Activate one conditional W/L guard per dynamic PET sequence browser."""
  if browser is None or sequence is None or proxy is None:
    return
  key = browser.GetID()
  previous = _dpetPETWindowGuards.get(key)
  if previous is not None:
    previous.close()
  _dpetPETWindowGuards[key] = _dpetPETSequenceWindowGuard(
    browser, sequence, proxy)


class dPETImporterPluginClass(DICOMPlugin):
  """ dPETImporter specific parser
  """

  def __init__(self,epsilon=0.01):
    super().__init__()
    self.loadType = "dPETImporter"

    self.tags['seriesInstanceUID'] = "0020,000E"
    self.tags['seriesDescription'] = "0008,103E"
    self.tags['instanceUID'] = "0008,0018"
    self.tags['position'] = "0020,0032"
    self.tags['orientation'] = "0020,0037"
    self.tags['studyDescription'] = "0008,1030"
    self.tags['seriesNumber'] = "0020,0011"
    self.tags['instanceNumber'] = "0020,0013"
    # self.tags['repetitionTime'] = "0018,0080"
    self.tags['modality'] = "0008,0060"
    self.tags['studyInstanceUID'] = "0020,000D"
    self.tags['frameOfReferenceUID'] = "0020,0052"
    self.tags['sopClassUID'] = "0008,0016"
    self.tags['seriesType'] = "0054,1000"
    self.tags['acquisitionDate'] = "0008,0022"
    self.tags['acquisitionDateTime'] = "0008,002A"
    self.tags['actualFrameDuration'] = "0018,1242"
    self.tags['units'] = "0054,1001"
    self.tags['suvType'] = "0054,1006"
    self.tags['decayCorrection'] = "0054,1102"
    self.tags['correctedImage'] = "0028,0051"

    # tags used to identify multivolumes
    self.multiVolumeTags = {
      'AcquisitionTime': "0008,0032",
      'SeriesTime': "0008,0031",
      'ContentTime': "0008,0033",
      'TemporalPositionIdentifier': "0020,0100",
      'TriggerTime': "0018,1060",
    }
    for k,v in self.multiVolumeTags.items():
      self.tags[k] = v

    self.multiVolumeTagsUnits = {
      'AcquisitionTime': "ms",
      'SeriesTime': "ms",
      'ContentTime': "ms",
      'TemporalPositionIdentifier': "count",
      'TriggerTime': "ms",
    }
    self.epsilon = epsilon
    self.detailedLogging = False

  @staticmethod
  def settingsPanelEntry(panel, parent):
    formLayout = qt.QFormLayout(parent)
    enabledCheckBox = qt.QCheckBox()
    enabledCheckBox.toolTip = "If enabled, dynamic PET series will be offered for loading as a volume sequence."
    enabledCheckBox.checked = True
    formLayout.addRow("Enable Dynamic PET importer:", enabledCheckBox)
    panel.registerProperty(
      "DICOM/dPETImporterEnabled",
      enabledCheckBox,
      "checked",
      str(qt.SIGNAL("toggled(bool)"))
    )
    staticCheckBox = qt.QCheckBox()
    staticCheckBox.toolTip = (
      "If enabled, STATIC and WHOLE BODY PET series are also offered as an "
      "optional metadata-preserving loadable for SlicerDynamicPET kinetic analysis. "
      "This alternative remains unchecked by default so specialized PET/SUV "
      "loaders can remain the preferred choice.")
    staticCheckBox.checked = True
    formLayout.addRow("Offer Static PET kinetic loader:", staticCheckBox)
    panel.registerProperty(
      "DICOM/dPETImporterStaticEnabled",
      staticCheckBox,
      "checked",
      str(qt.SIGNAL("toggled(bool)"))
    )
    suvCheckBox = qt.QCheckBox()
    suvCheckBox.toolTip = "If enabled, each frame will be converted to SUVbw (if metadata is sufficient)."
    suvCheckBox.checked = True
    formLayout.addRow("Convert to SUVbw on load:", suvCheckBox)
    panel.registerProperty(
      "DICOM/dPETImporterSUVEnabled",
      suvCheckBox,
      "checked",
      str(qt.SIGNAL("toggled(bool)"))
    )

  def _missingSUVKeys(self, mvNode):
    required = [
      'RadionuclideStartDateTime',
      'RadionuclideHalfLife',
      'RadionuclideTotalDose',
      'PatientWeight',
      'DecayCorrection',
      'CorrectedImage',
      'dPET.FirstFrameAcquisitionDateTime',
    ]
    missing = []
    for k in required:
      v = mvNode.GetAttribute(k) if mvNode else None
      if v is None or v == "":
        missing.append(k)
    return missing

  def _isDecayCorrectedToFirstFrameStart(self, mvNode):
    """
    Enforce our convention:
      - DICOM DecayCorrection must be START
      - first frame acquisition datetime must be known
    """
    if mvNode is None:
      return False
    decayCorrection = (mvNode.GetAttribute('DecayCorrection') or '').upper()
    firstFrameDt = mvNode.GetAttribute('dPET.FirstFrameAcquisitionDateTime')
    return (decayCorrection == "START") and bool(firstFrameDt)

  def _getRadiopharmNested(self, dicomFilePath):
    """
    Read nested Radiopharmaceutical Information Sequence (0054,0016)
    using pydicom, return dict with keys:
      RadionuclideHalfLife, RadionuclideTotalDose,
      RadiopharmaceuticalStartDateTime, RadiopharmaceuticalStartTime
    Returns {} if not found.
    """
    try:
      import pydicom
      ds = pydicom.dcmread(dicomFilePath, stop_before_pixels=True, force=True)
      de = ds.get((0x0054, 0x0016), None)  # RadiopharmaceuticalInformationSequence
      if de is None:
        return {}

      seq = getattr(de, "value", de)
      if seq is None:
        return {}
      if not hasattr(seq, "__len__"):
        # Unexpected type
        logging.error(f"[dPET] Radiopharm sequence has unexpected type: {type(seq)}")
        return {}
      if len(seq) == 0:
        return {}

      item = seq[0]

      out = {}
      def get_item_str(tag):
        if tag in item:
          v = item[tag].value
          return "" if v is None else str(v)
        return ""
      out["RadionuclideHalfLife"] = get_item_str((0x0018, 0x1075))
      out["RadionuclideTotalDose"] = get_item_str((0x0018, 0x1074))
      out["RadiopharmaceuticalStartDateTime"] = get_item_str((0x0018, 0x1078))
      out["RadiopharmaceuticalStartTime"] = get_item_str((0x0018, 0x1072))
      out["RadiopharmaceuticalName"] = get_item_str((0x0018, 0x0031))

      def get_code_sequence(tag):
        if tag not in item:
          return ""
        seq_value = getattr(item[tag], "value", None)
        if not seq_value or len(seq_value) == 0:
          return ""
        code = seq_value[0]
        cv = str(getattr(code, "CodeValue", "") or "")
        csd = str(getattr(code, "CodingSchemeDesignator", "") or "")
        cm = str(getattr(code, "CodeMeaning", "") or "")
        if not (cv or csd or cm):
          return ""
        return "|".join((cv, csd, cm))

      out["RadionuclideCode"] = get_code_sequence((0x0054, 0x0300))
      out["RadiopharmaceuticalCode"] = get_code_sequence((0x0054, 0x0304))

      out = {k: v for k, v in out.items() if v not in ("", None)}
      return out
    except Exception as e:
      logging.error(f"[dPET] pydicom nested radiopharm read failed: {e}")
      import traceback
      traceback.print_exc()
      return {}


  def _resolveRadiopharmaceuticalStartDateTime(self, mvNode, firstFrameDt):
    """
    Resolve the best available administration datetime while preserving source.

    Priority:
      1. DICOM Radiopharmaceutical Start DateTime (0018,1078)
      2. DICOM Radiopharmaceutical Start Time (0018,1072) combined with
         the first-frame calendar date (with +/- one-day rollover candidates)

    The resolved datetime is metadata only.  DynamicPET applies a separate,
    conservative rule before using it as a kinetic time offset.
    """
    exact = mvNode.GetAttribute('RadiopharmaceuticalStartDateTime') or ''
    exactDt = _parse_dicom_datetime_value(exact)
    if exactDt is not None:
      return exactDt, 'DICOM.0018,1078'

    startTime = mvNode.GetAttribute('RadiopharmaceuticalStartTime') or ''
    if not startTime or firstFrameDt is None:
      return None, 'Unavailable'

    # Parse TM using the first-frame date, then consider midnight rollover.
    baseDate = firstFrameDt.strftime('%Y%m%d')
    candidate = _parse_dicom_datetime(baseDate, startTime)
    if candidate is None:
      return None, 'Unavailable'

    candidates = [candidate - timedelta(days=1), candidate, candidate + timedelta(days=1)]
    # Choose the closest calendar interpretation to the first PET frame.
    resolved = min(candidates, key=lambda dt: abs((firstFrameDt - dt).total_seconds()))
    return resolved, 'DICOM.0018,1072+frame-date'


  @staticmethod
  def _formatDicomDateTime(dt):
    if dt is None:
      return ""
    return dt.strftime("%Y%m%d%H%M%S.%f").rstrip('0').rstrip('.')

  def _fileValue(self, filePath, tag):
    try:
      value = slicer.dicomDatabase.fileValue(filePath, tag)
      return "" if value in (None, "") else str(value)
    except Exception:
      return ""

  def _isPETImageFile(self, filePath):
    """Return True only for DICOM objects whose Modality (0008,0060) is PT.

    dPETImporter is an image importer. Derived DICOM objects that can coexist
    with PET series in the same study (for example Real World Value Mapping,
    RTSTRUCT, Parametric Map, CT, etc.) must be ignored before any dynamic-PET
    grouping or PET-specific metadata parsing is attempted.
    """
    return self._fileValue(filePath, self.tags['modality']).strip().upper() == 'PT'

  def _petImageFiles(self, files):
    """Silently remove non-PET DICOM objects from an examine/load candidate."""
    return [filePath for filePath in (files or []) if self._isPETImageFile(filePath)]

  def _acquisitionDateTimeForFile(self, filePath):
    """Return (datetime, source) for the real acquisition start of one PET image."""
    acqDateTime = self._fileValue(filePath, "0008,002A")
    parsed = _parse_dicom_datetime_value(acqDateTime)
    if parsed is not None:
      return parsed, "DICOM.0008,002A"

    acqDate = self._fileValue(filePath, "0008,0022")
    acqTime = self._fileValue(filePath, "0008,0032")
    parsed = _parse_dicom_datetime(acqDate, acqTime)
    if parsed is not None:
      return parsed, "DICOM.0008,0022+0008,0032"

    # Last-resort series clock. Keep the source explicit because this is less
    # specific than Acquisition Date/Time and should be treated conservatively.
    seriesDate = self._fileValue(filePath, "0008,0021") or self._fileValue(filePath, "0008,0020")
    seriesTime = self._fileValue(filePath, "0008,0031")
    parsed = _parse_dicom_datetime(seriesDate, seriesTime)
    if parsed is not None:
      return parsed, "DICOM.series-date-time-fallback"
    return None, "Unavailable"

  def _actualFrameDurationSecForFile(self, filePath):
    """Return acquisition duration in seconds, or None if unavailable."""
    value = self._fileValue(filePath, "0018,1242")  # Actual Frame Duration, ms
    if value:
      try:
        duration = float(str(value).split('\\')[0]) / 1000.0
        if duration > 0.0:
          return duration, "DICOM.0018,1242"
      except Exception:
        pass

    # Existing dPET convention used by some test/vendor data: seconds.
    value = self._fileValue(filePath, "0067,1004")
    if value:
      try:
        duration = float(str(value).split('\\')[0])
        if duration > 0.0:
          return duration, "DICOM.0067,1004"
      except Exception:
        pass
    return None, "Unavailable"

  def _seriesTypeValue1(self, filePath):
    value = self._fileValue(filePath, "0054,1000").upper().strip()
    if not value:
      return ""
    return value.split('\\')[0].strip()

  def _resolveRadiopharmaceuticalStartDateTimeFromValues(
      self, exactValue, startTimeValue, firstAcquisitionDt):
    exactDt = _parse_dicom_datetime_value(exactValue or '')
    if exactDt is not None:
      return exactDt, 'DICOM.0018,1078'

    startTime = str(startTimeValue or '').strip()
    if not startTime or firstAcquisitionDt is None:
      return None, 'Unavailable'

    candidate = _parse_dicom_datetime(firstAcquisitionDt.strftime('%Y%m%d'), startTime)
    if candidate is None:
      return None, 'Unavailable'

    candidates = [candidate - timedelta(days=1), candidate, candidate + timedelta(days=1)]
    resolved = min(candidates, key=lambda dt: abs((firstAcquisitionDt - dt).total_seconds()))
    return resolved, 'DICOM.0018,1072+acquisition-date'

  def _commonPETMetadata(self, firstFile, firstAcquisitionDt):
    """Collect compact quantitative/radiopharmaceutical metadata at import time."""
    directTags = {
      'PatientID': '0010,0020',
      'PatientWeight': '0010,1030',
      'Units': '0054,1001',
      'SUVType': '0054,1006',
      'DecayCorrection': '0054,1102',
      'CorrectedImage': '0028,0051',
      'DecayFactor': '0054,1321',
      'FrameReferenceTime': '0054,1300',
      'RadionuclideHalfLife': '0018,1075',
      'RadionuclideTotalDose': '0018,1074',
      'RadiopharmaceuticalStartDateTime': '0018,1078',
      'RadiopharmaceuticalStartTime': '0018,1072',
      'RadiopharmaceuticalName': '0018,0031',
    }
    out = {}
    for name, tag in directTags.items():
      value = self._fileValue(firstFile, tag)
      if value:
        out[name] = value

    nested = self._getRadiopharmNested(firstFile)
    for key, value in nested.items():
      if value and not out.get(key):
        out[key] = str(value)

    startDt, startSource = self._resolveRadiopharmaceuticalStartDateTimeFromValues(
      out.get('RadiopharmaceuticalStartDateTime', ''),
      out.get('RadiopharmaceuticalStartTime', ''),
      firstAcquisitionDt)
    if startDt is not None:
      out['RadionuclideStartDateTime'] = self._formatDicomDateTime(startDt)
      out['InjectionDateTimeSource'] = startSource
      if firstAcquisitionDt is not None:
        out['InjectionToAcquisitionOffsetSec'] = float((firstAcquisitionDt - startDt).total_seconds())
    else:
      out['InjectionDateTimeSource'] = 'Unavailable'
    return out

  def _buildStaticKineticMetadata(self, files):
    """Build self-contained timing/provenance metadata for STATIC/WHOLE BODY PET."""
    if not files:
      return None

    firstFile = files[0]
    seriesType = self._seriesTypeValue1(firstFile)
    if seriesType not in ('STATIC', 'WHOLE BODY'):
      return None

    records = []
    startDts = []
    endDts = []
    timingKeys = set()
    startKeys = set()
    durationKeys = set()
    timingComplete = True
    acquisitionSources = set()
    durationSources = set()

    for filePath in files:
      startDt, startSource = self._acquisitionDateTimeForFile(filePath)
      durationSec, durationSource = self._actualFrameDurationSecForFile(filePath)
      acquisitionSources.add(startSource)
      durationSources.add(durationSource)

      if startDt is None or durationSec is None:
        timingComplete = False
      if startDt is not None:
        startDts.append(startDt)
        startKeys.add(self._formatDicomDateTime(startDt))
      if durationSec is not None:
        durationKeys.add(round(float(durationSec), 6))
      endDt = None
      if startDt is not None and durationSec is not None:
        endDt = startDt + timedelta(seconds=float(durationSec))
        endDts.append(endDt)
        timingKeys.add((self._formatDicomDateTime(startDt), round(float(durationSec), 6)))

      record = {
        'sopInstanceUID': self._fileValue(filePath, '0008,0018'),
        'sopClassUID': self._fileValue(filePath, '0008,0016'),
        'instanceNumber': self._fileValue(filePath, '0020,0013'),
        'acquisitionStartDateTime': self._formatDicomDateTime(startDt),
        'acquisitionStartSource': startSource,
        'durationSec': durationSec,
        'durationSource': durationSource,
        'acquisitionEndDateTime': self._formatDicomDateTime(endDt),
      }
      position = self._parseDICOMVector(self._fileValue(filePath, '0020,0032'), 3)
      orientation = self._parseDICOMVector(self._fileValue(filePath, '0020,0037'), 6)
      if position is not None:
        record['imagePositionPatient'] = position
      if orientation is not None:
        record['imageOrientationPatient'] = orientation
      records.append(record)

    earliestStart = min(startDts) if startDts else None
    latestEnd = max(endDts) if endDts else None
    spatiallyVaryingTiming = (len(startKeys) > 1 or len(durationKeys) > 1 or len(timingKeys) > 1)
    timingMode = 'INCOMPLETE'
    if timingComplete:
      timingMode = 'SPATIAL' if spatiallyVaryingTiming else 'UNIFORM'

    common = self._commonPETMetadata(firstFile, earliestStart)
    metadata = {
      'schemaVersion': 1,
      'metadataSource': 'dPETImporter',
      'acquisitionKind': 'STATIC',
      'seriesType': seriesType,
      'wholeBody': bool(seriesType == 'WHOLE BODY' or spatiallyVaryingTiming),
      'timingMode': timingMode,
      'timingComplete': bool(timingComplete),
      'spatiallyVaryingTiming': bool(spatiallyVaryingTiming),
      'acquisitionStartDateTime': self._formatDicomDateTime(earliestStart),
      'acquisitionEndDateTime': self._formatDicomDateTime(latestEnd),
      'timingSources': sorted(acquisitionSources),
      'durationSources': sorted(durationSources),
      'studyInstanceUID': self._fileValue(firstFile, '0020,000D'),
      'seriesInstanceUID': self._fileValue(firstFile, '0020,000E'),
      'frameOfReferenceUID': self._fileValue(firstFile, '0020,0052'),
      'common': common,
      'spatialTiming': records,
    }
    return metadata

  def _buildDynamicKineticMetadata(self, files, nFrames, filesPerFrame, mvNode):
    """Build the same metadata contract for a dPETImporter dynamic sequence."""
    if not files or nFrames < 1 or filesPerFrame < 1:
      return None
    firstFile = files[0]
    frames = []
    startDts = []
    endDts = []
    timingComplete = True
    for frameIndex in range(nFrames):
      filePath = files[frameIndex * filesPerFrame]
      startDt, startSource = self._acquisitionDateTimeForFile(filePath)
      durationSec, durationSource = self._actualFrameDurationSecForFile(filePath)
      if startDt is None or durationSec is None:
        timingComplete = False
      if startDt is not None:
        startDts.append(startDt)
      endDt = None
      if startDt is not None and durationSec is not None:
        endDt = startDt + timedelta(seconds=float(durationSec))
        endDts.append(endDt)
      frames.append({
        'index': frameIndex,
        'acquisitionStartDateTime': self._formatDicomDateTime(startDt),
        'acquisitionStartSource': startSource,
        'durationSec': durationSec,
        'durationSource': durationSource,
        'acquisitionEndDateTime': self._formatDicomDateTime(endDt),
      })

    earliestStart = min(startDts) if startDts else None
    latestEnd = max(endDts) if endDts else None
    common = self._commonPETMetadata(firstFile, earliestStart)

    # Prefer already normalized values from the existing dynamic parser.
    for attr, key in (
      ('RadionuclideStartDateTime', 'RadionuclideStartDateTime'),
      ('dPET.InjectionDateTimeSource', 'InjectionDateTimeSource'),
      ('dPET.InjectionToAcquisitionOffsetSec', 'InjectionToAcquisitionOffsetSec'),
    ):
      value = mvNode.GetAttribute(attr) if mvNode else None
      if value not in (None, ''):
        if key == 'InjectionToAcquisitionOffsetSec':
          try:
            value = float(value)
          except Exception:
            pass
        common[key] = value

    return {
      'schemaVersion': 1,
      'metadataSource': 'dPETImporter',
      'acquisitionKind': 'DYNAMIC',
      'seriesType': self._seriesTypeValue1(firstFile) or 'DYNAMIC',
      'wholeBody': False,
      'timingMode': 'FRAMES',
      'timingComplete': bool(timingComplete),
      'spatiallyVaryingTiming': False,
      'acquisitionStartDateTime': self._formatDicomDateTime(earliestStart),
      'acquisitionEndDateTime': self._formatDicomDateTime(latestEnd),
      'studyInstanceUID': self._fileValue(firstFile, '0020,000D'),
      'seriesInstanceUID': self._fileValue(firstFile, '0020,000E'),
      'frameOfReferenceUID': self._fileValue(firstFile, '0020,0052'),
      'common': common,
      'frames': frames,
    }

  def _applyKineticMetadata(self, node, metadata):
    if node is None or not metadata:
      return
    compact = json.dumps(metadata, separators=(',', ':'), ensure_ascii=False)
    node.SetAttribute('dPET.KineticMetadataSchemaVersion', str(metadata.get('schemaVersion', 1)))
    node.SetAttribute('dPET.KineticMetadata', compact)
    node.SetAttribute('dPET.AcquisitionKind', str(metadata.get('acquisitionKind', '')))
    node.SetAttribute('dPET.SeriesType', str(metadata.get('seriesType', '')))
    node.SetAttribute('dPET.AcquisitionTimingMode', str(metadata.get('timingMode', '')))
    node.SetAttribute('dPET.AcquisitionTimingComplete', '1' if metadata.get('timingComplete') else '0')
    node.SetAttribute('dPET.SpatiallyVaryingTiming', '1' if metadata.get('spatiallyVaryingTiming') else '0')
    node.SetAttribute('dPET.WholeBody', '1' if metadata.get('wholeBody') else '0')

    startText = metadata.get('acquisitionStartDateTime') or ''
    endText = metadata.get('acquisitionEndDateTime') or ''
    if startText:
      node.SetAttribute('dPET.AcquisitionStartDateTime', str(startText))
      node.SetAttribute('dPET.FirstFrameAcquisitionDateTime', str(startText))
    if endText:
      node.SetAttribute('dPET.AcquisitionEndDateTime', str(endText))

    for key, attr in (
      ('studyInstanceUID', 'dPET.DICOM.StudyInstanceUID'),
      ('seriesInstanceUID', 'dPET.DICOM.SeriesInstanceUID'),
      ('frameOfReferenceUID', 'dPET.DICOM.FrameOfReferenceUID'),
    ):
      value = metadata.get(key)
      if value:
        node.SetAttribute(attr, str(value))

    if metadata.get('acquisitionKind') == 'STATIC':
      node.SetAttribute(
        'dPET.Static.SpatialTiming',
        json.dumps(metadata.get('spatialTiming', []), separators=(',', ':'), ensure_ascii=False))

    common = metadata.get('common') or {}
    commonAttributeMap = {
      'PatientID': 'dPET.PatientID',
      'PatientWeight': 'PatientWeight',
      'Units': 'Units',
      'SUVType': 'SUVType',
      'DecayCorrection': 'DecayCorrection',
      'CorrectedImage': 'CorrectedImage',
      'DecayFactor': 'DecayFactor',
      'FrameReferenceTime': 'FrameReferenceTime',
      'RadionuclideHalfLife': 'RadionuclideHalfLife',
      'RadionuclideTotalDose': 'RadionuclideTotalDose',
      'RadiopharmaceuticalStartDateTime': 'RadiopharmaceuticalStartDateTime',
      'RadiopharmaceuticalStartTime': 'RadiopharmaceuticalStartTime',
      'RadiopharmaceuticalName': 'RadiopharmaceuticalName',
      'RadionuclideCode': 'dPET.RadionuclideCode',
      'RadiopharmaceuticalCode': 'dPET.RadiopharmaceuticalCode',
      'RadionuclideStartDateTime': 'RadionuclideStartDateTime',
      'InjectionDateTimeSource': 'dPET.InjectionDateTimeSource',
      'InjectionToAcquisitionOffsetSec': 'dPET.InjectionToAcquisitionOffsetSec',
    }
    for key, attr in commonAttributeMap.items():
      value = common.get(key)
      if value not in (None, ''):
        node.SetAttribute(attr, str(value))

  def examineStaticFiles(self, files):
    """Offer one optional metadata-preserving scalar loadable per static PET series."""
    if not settingsValue('DICOM/dPETImporterStaticEnabled', True, converter=toBool):
      return []
    if not files:
      return []

    seriesLists = {}
    for filePath in files:
      if self._fileValue(filePath, self.tags['modality']).upper() != 'PT':
        continue
      sid = self._fileValue(filePath, self.tags['seriesInstanceUID']) or 'Unknown'
      seriesLists.setdefault(sid, []).append(filePath)

    scalarPluginClass = slicer.modules.dicomPlugins.get('DICOMScalarVolumePlugin')
    if scalarPluginClass is None:
      return []
    scalarPlugin = scalarPluginClass()

    loadables = []
    for sid, seriesFiles in seriesLists.items():
      seriesType = self._seriesTypeValue1(seriesFiles[0])
      if seriesType not in ('STATIC', 'WHOLE BODY'):
        continue

      scalarLoadables = scalarPlugin.examine([seriesFiles]) or []
      scalarLoadables = [item for item in scalarLoadables if getattr(item, 'files', None)]
      if not scalarLoadables:
        continue
      # Prefer the interpretation that retains the complete series.
      scalarLoadable = max(scalarLoadables, key=lambda item: len(item.files))
      metadata = self._buildStaticKineticMetadata(scalarLoadable.files)
      if not metadata:
        continue

      loadable = DICOMLoadable()
      loadable.files = list(scalarLoadable.files)
      baseName = scalarLoadable.name or self._fileValue(seriesFiles[0], self.tags['seriesDescription']) or 'Static PET'
      loadable.name = f'{baseName} [dPET kinetic metadata]'
      timingMode = metadata.get('timingMode', 'INCOMPLETE')
      if timingMode == 'SPATIAL':
        timingText = 'spatially varying acquisition timing preserved'
      elif timingMode == 'UNIFORM':
        timingText = 'uniform acquisition timing preserved'
      else:
        timingText = 'acquisition timing incomplete'
      loadable.tooltip = (
        'Optional SlicerDynamicPET static PET loader; keeps kinetic timing, radiopharmaceutical, '
        f'quantitative, and DICOM identity metadata in MRML ({timingText}). '
        'Left unchecked by default so dedicated PET/SUV importers remain preferred.')
      loadable.selected = False
      loadable.confidence = 0.40
      loadable.dPETLoadMode = 'static'
      loadable.dPETScalarLoadable = scalarLoadable
      loadable.dPETKineticMetadata = metadata
      loadables.append(loadable)
    return loadables


  def _resolveEnhancedQuantification(self, ds, frameGroup, sharedGroup):
    """Return (units, slope, intercept, mapped_min, mapped_max)."""
    def first(group, name):
      items = getattr(group, name, None) if group is not None else None
      return items[0] if items else None

    rwvm = (getattr(frameGroup, 'RealWorldValueMappingSequence', None)
            or getattr(sharedGroup, 'RealWorldValueMappingSequence', None)
            or getattr(ds, 'RealWorldValueMappingSequence', None))
    if rwvm:
      matches = []
      for item in rwvm:
        code = first(item, 'MeasurementUnitsCodeSequence')
        value = str(getattr(code, 'CodeValue', '')).upper().replace(' ', '')
        scheme = str(getattr(code, 'CodingSchemeDesignator', '')).upper()
        units = ('BQML' if scheme == 'UCUM' and value in ('BQ/ML', 'BQ.ML-1')
                 else 'GML' if scheme == 'UCUM' and value in ('G/ML', 'G.ML-1')
                 else None)
        if units:
          if getattr(item, 'RealWorldValueLUTData', None) is not None:
            raise ValueError('Enhanced PET non-linear real-world mapping is unsupported')
          slope = float(item.RealWorldValueSlope)
          intercept = float(item.RealWorldValueIntercept)
          low = getattr(item, 'DoubleFloatRealWorldValueFirstValueMapped', None)
          high = getattr(item, 'DoubleFloatRealWorldValueLastValueMapped', None)
          if low is None:
            low = getattr(item, 'RealWorldValueFirstValueMapped', None)
          if high is None:
            high = getattr(item, 'RealWorldValueLastValueMapped', None)
          matches.append((units, slope, intercept,
                          float(low) if low is not None else None,
                          float(high) if high is not None else None))
      # Prefer activity concentration when both activity and SUV maps exist.
      matches.sort(key=lambda x: x[0] != 'BQML')
      if not matches:
        raise ValueError('Enhanced PET mapping lacks supported Bq/mL or g/mL UCUM units')
      if sum(m[0] == matches[0][0] for m in matches) != 1:
        raise ValueError('Ambiguous Enhanced PET real-world mappings for the same units')
      mapping = matches[0]
    else:
      # Fall back to rescale only if its physical units are explicitly known.
      frameType = first(frameGroup, 'PETFrameTypeSequence') or first(
        sharedGroup, 'PETFrameTypeSequence')
      units = str(getattr(ds, 'Units', None) or getattr(frameType, 'Units', '')).upper()
      pvt = first(frameGroup, 'PixelValueTransformationSequence') or first(
        sharedGroup, 'PixelValueTransformationSequence')
      if units not in ('BQML', 'GML') or pvt is None:
        raise ValueError('Enhanced PET has no usable quantitative mapping with units')
      mapping = (units, float(pvt.RescaleSlope),
                 float(pvt.RescaleIntercept), None, None)

    import math
    if (not math.isfinite(mapping[1]) or not math.isfinite(mapping[2])
        or mapping[1] == 0):
      raise ValueError('Invalid Enhanced PET pixel-value transformation')
    low, high = mapping[3:]
    if ((low is not None and not math.isfinite(low))
        or (high is not None and not math.isfinite(high))
        or (low is not None and high is not None and low > high)):
      raise ValueError('Invalid Enhanced PET real-world mapped pixel range')
    frameType = first(frameGroup, 'PETFrameTypeSequence') or first(
      sharedGroup, 'PETFrameTypeSequence')
    suvType = str(getattr(ds, 'SUVType', None) or getattr(frameType, 'SUVType', None) or 'BW').upper()
    if mapping[0] == 'GML' and suvType not in ('BW', 'SUVBW'):
      raise ValueError('Only SUVbw is supported (GML with BW normalization)')
    return mapping

  def _readEnhancedFrameMetadata(self, filePath):
    """Validate one Enhanced PET file without reading its pixel data."""
    import numpy as np
    import pydicom
    ds = pydicom.dcmread(filePath, stop_before_pixels=True)
    if str(getattr(ds, 'SOPClassUID', '')) != '1.2.840.10008.5.1.4.1.1.130':
      return None

    shared = (getattr(ds, 'SharedFunctionalGroupsSequence', None) or [None])[0]
    frameGroups = getattr(ds, 'PerFrameFunctionalGroupsSequence', None)
    if not frameGroups or len(frameGroups) != int(getattr(ds, 'NumberOfFrames', 0)):
      raise ValueError('Enhanced PET has missing/mismatched per-frame functional groups')
    if int(getattr(ds, 'SamplesPerPixel', 0)) != 1:
      raise ValueError('Only single-channel Enhanced PET images are supported')
    if str(getattr(ds, 'PhotometricInterpretation', '')) != 'MONOCHROME2':
      raise ValueError('Only MONOCHROME2 Enhanced PET images are supported')
    rows, cols = int(ds.Rows), int(ds.Columns)
    if rows <= 0 or cols <= 0:
      raise ValueError('Invalid Enhanced PET image dimensions')
    paddingValue = getattr(ds, 'PixelPaddingValue', None)
    paddingLimit = getattr(ds, 'PixelPaddingRangeLimit', None)
    if paddingLimit is not None and paddingValue is None:
      raise ValueError('Enhanced PET PixelPaddingRangeLimit has no PixelPaddingValue')

    def first(group, name):
      items = getattr(group, name, None) if group is not None else None
      return items[0] if items else None

    frames = {}
    sourceToFrame = []
    mappings = []
    firstSeen = []
    for index, fg in enumerate(frameGroups):
      fc = first(fg, 'FrameContentSequence')
      if fc is None or getattr(fc, 'TemporalPositionIndex', None) is None:
        raise ValueError('FrameContentSequence.TemporalPositionIndex missing')
      t = int(fc.TemporalPositionIndex)
      stack = str(getattr(fc, 'StackID', ''))
      sliceIndex = getattr(fc, 'InStackPositionNumber', None)
      if not stack or sliceIndex is None or t < 1 or int(sliceIndex) < 1:
        raise ValueError('Enhanced PET stack or in-stack index missing/invalid')
      if frames and stack != frames[firstSeen[0]]['stack']:
        raise ValueError('Multiple Enhanced PET stacks are not supported')
      pm = first(fg, 'PixelMeasuresSequence') or first(shared, 'PixelMeasuresSequence')
      po = first(fg, 'PlaneOrientationSequence') or first(shared, 'PlaneOrientationSequence')
      pp = first(fg, 'PlanePositionSequence') or first(shared, 'PlanePositionSequence')
      if pm is None or po is None or pp is None:
        raise ValueError('Enhanced PET lacks spatial functional groups')
      spacing = np.asarray(pm.PixelSpacing, dtype=float)
      orient = np.asarray(po.ImageOrientationPatient, dtype=float)
      position = np.asarray(pp.ImagePositionPatient, dtype=float)
      if spacing.shape != (2,) or orient.shape != (6,) or position.shape != (3,):
        raise ValueError('Invalid Enhanced PET spatial vectors')
      if not (np.isfinite(spacing).all() and np.isfinite(orient).all()
              and np.isfinite(position).all() and np.all(spacing > 0)):
        raise ValueError('Non-finite Enhanced PET geometry')
      start = _parse_dicom_datetime_value(getattr(fc, 'FrameAcquisitionDateTime', None))
      if start is None:
        raise ValueError('Enhanced PET frame acquisition datetime missing/invalid')
      duration = getattr(fc, 'FrameAcquisitionDuration', None)
      duration = float(duration) / 1000.0 if duration is not None else None
      if duration is not None and (not np.isfinite(duration) or duration <= 0):
        raise ValueError('Invalid Enhanced PET frame acquisition duration')
      mapping = self._resolveEnhancedQuantification(ds, fg, shared)
      mappings.append(mapping)
      if t not in frames:
        frames[t] = {'stack': stack, 'start': start, 'duration': duration,
                     'indices': [], 'slices': {}, 'spacing': spacing, 'orientation': orient,
                     'unit': mapping[0]}
        firstSeen.append(t)
      record = frames[t]
      k = int(sliceIndex)
      if k in record['slices']:
        raise ValueError(f'Duplicate in-stack slice {k} in temporal position {t}')
      if (record['stack'] != stack or record['start'] != start
          or record['duration'] != duration or record['unit'] != mapping[0]
          or not np.allclose(record['spacing'], spacing, atol=1e-5)
          or not np.allclose(record['orientation'], orient, atol=1e-5)):
        raise ValueError(f'Inconsistent geometry, timing or units within temporal position {t}')
      record['indices'].append(index)
      record['slices'][k] = (index, position)
      sourceToFrame.append(t)

    if len(frames) < 2:
      return None  # Existing single-volume importers remain responsible.
    temporalOrdinals = {value: i for i, value in enumerate(firstSeen)}
    sourceToFrame = [temporalOrdinals[t] for t in sourceToFrame]
    if firstSeen != sorted(firstSeen):
      raise ValueError('Enhanced PET temporal positions are not in ascending order')
    if any(sourceToFrame[i] < sourceToFrame[i-1] for i in range(1, len(sourceToFrame))):
      raise ValueError('Interleaved Enhanced PET time points need a different decoder path')

    reference = frames[firstSeen[0]]
    ref_keys = set(reference['slices'])
    if len(ref_keys) < 2:
      raise ValueError('Enhanced PET volume requires at least two spatial slices')
    orient = reference['orientation']
    axisI, axisJ = orient[:3], orient[3:]
    normal = np.cross(axisI, axisJ)
    if (not np.isclose(np.linalg.norm(axisI), 1, atol=1e-4)
        or not np.isclose(np.linalg.norm(axisJ), 1, atol=1e-4)
        or not np.isclose(np.linalg.norm(normal), 1, atol=1e-4)):
      raise ValueError('Invalid Enhanced PET orientation vectors')
    ref_sorted = sorted(ref_keys, key=lambda k: np.dot(reference['slices'][k][1], normal))
    ref_positions = np.asarray([reference['slices'][k][1] for k in ref_sorted])
    sliceVector = (ref_positions[-1] - ref_positions[0]) / (len(ref_sorted) - 1)
    if np.dot(sliceVector, normal) <= 0 or np.linalg.norm(sliceVector) < 1e-5:
      raise ValueError('Invalid Enhanced PET slice progression')
    for i, pos in enumerate(ref_positions):
      if not np.allclose(pos, ref_positions[0] + i * sliceVector, atol=0.05):
        raise ValueError('Enhanced PET reference slices are not evenly spaced')
    units = reference['unit']
    starts = []
    durations = []
    sourceToK = [None] * len(frameGroups)
    for fi, t in enumerate(firstSeen):
      frame = frames[t]
      if (set(frame['slices']) != ref_keys or frame['unit'] != units
          or not np.allclose(frame['spacing'], reference['spacing'], atol=1e-5)
          or not np.allclose(frame['orientation'], orient, atol=1e-5)):
        raise ValueError(f'Enhanced PET temporal position {t} has incompatible slices')
      for kIndex, key in enumerate(ref_sorted):
        src, pos = frame['slices'][key]
        if not np.allclose(pos, ref_positions[kIndex], atol=0.05):
          raise ValueError(f'Enhanced PET temporal position {t} changes spatial geometry')
        sourceToK[src] = kIndex
      starts.append(frame['start'])
      durations.append(frame['duration'])
    if any(starts[i] <= starts[i-1] for i in range(1, len(starts))):
      raise ValueError('Enhanced PET acquisition times are not strictly increasing')

    # Build IJK-to-RAS using DICOM LPS row, column and measured K vectors.
    spacing = reference['spacing']
    lpsMatrix = np.eye(4)
    lpsMatrix[:3, 0] = axisI * spacing[1]
    lpsMatrix[:3, 1] = axisJ * spacing[0]
    lpsMatrix[:3, 2] = sliceVector
    lpsMatrix[:3, 3] = ref_positions[0]
    rasMatrix = np.diag([-1., -1., 1., 1.]) @ lpsMatrix

    frameType = first(frameGroups[0], 'PETFrameTypeSequence') or first(
      shared, 'PETFrameTypeSequence')
    return {
      'nFrames': len(firstSeen), 'nSlices': len(ref_sorted),
      'rows': rows, 'cols': cols, 'units': units,
      'times': starts, 'durations': durations, 'mappings': mappings,
      'sourceToK': sourceToK, 'sourceToFrame': sourceToFrame,
      'sourceFrameNumbers': [[frame['slices'][k][0] + 1 for k in ref_sorted]
                             for frame in (frames[t] for t in firstSeen)],
      'rasMatrix': rasMatrix, 'sopInstanceUID': str(ds.SOPInstanceUID),
      'sopClassUID': str(ds.SOPClassUID),
      'studyInstanceUID': str(getattr(ds, 'StudyInstanceUID', '')),
      'seriesInstanceUID': str(getattr(ds, 'SeriesInstanceUID', '')),
      'frameOfReferenceUID': str(getattr(ds, 'FrameOfReferenceUID', '')),
      'paddingValue': int(paddingValue) if paddingValue is not None else None,
      'paddingLimit': int(paddingLimit) if paddingLimit is not None else None,
      'decayCorrection': str(getattr(ds, 'DecayCorrection', '')),
      'suvType': str(getattr(ds, 'SUVType', None) or getattr(frameType, 'SUVType', None) or 'BW'),
    }

  def _enhancedQuantitativeMetadata(self, filePath, firstFrameStart):
    """Translate Enhanced PET correction/dose metadata to the classic MRML contract."""
    import math
    import pydicom
    ds = pydicom.dcmread(filePath, stop_before_pixels=True)
    metadata = self._commonPETMetadata(filePath, firstFrameStart)
    extra = {}

    # Enhanced PET isotope dose is encoded in MBq; the classic helper expects Bq.
    agents = getattr(ds, 'RadiopharmaceuticalInformationSequence', None) or []
    doseBq = None
    if len(agents) == 1:
      try:
        doseMBq = float(agents[0].RadionuclideTotalDose)
        if math.isfinite(doseMBq) and 0 < doseMBq <= 10000:
          doseBq = doseMBq * 1e6
          metadata['RadionuclideTotalDose'] = str(doseBq)
          extra['dPET.Enhanced.DoseSourceUnit'] = 'MBq'
          extra['dPET.Enhanced.DoseSourceValue'] = str(doseMBq)
      except (AttributeError, TypeError, ValueError, OverflowError):
        pass
    if doseBq is None:
      metadata.pop('RadionuclideTotalDose', None)
      logging.warning('[dPET Enhanced] Dose unavailable, ambiguous, or outside expected MBq range')

    # Classic CorrectedImage is represented by individual Enhanced PET flags.
    correctionTags = (
      ('DecayCorrected', 'DECY'), ('AttenuationCorrected', 'ATTN'),
      ('ScatterCorrected', 'SCAT'), ('DeadTimeCorrected', 'DTIM'),
      ('RandomsCorrected', 'RAN'), ('CountLossNormalizationCorrected', 'CLN'),
      ('SensitivityCalibrated', 'DCAL'), ('PatientMotionCorrected', 'MOTN'),
    )
    corrected = [code for tag, code in correctionTags
                 if str(getattr(ds, tag, '')).upper() == 'YES']
    if corrected:
      metadata['CorrectedImage'] = '\\'.join(corrected)
      extra['dPET.Enhanced.CorrectedImageDerived'] = '1'
    for tag, _ in correctionTags:
      text = str(getattr(ds, tag, '') or '').upper()
      if text in ('YES', 'NO'):
        extra['dPET.Enhanced.' + tag] = text

    decay = str(getattr(ds, 'DecayCorrected', '')).upper()
    referenceText = str(getattr(ds, 'DecayCorrectionDateTime', '') or '')
    reference = _parse_dicom_datetime_value(referenceText)
    if referenceText:
      extra['dPET.DecayCorrectionDateTime'] = referenceText
    injected = _parse_dicom_datetime_value(
      metadata.get('RadionuclideStartDateTime') or
      metadata.get('RadiopharmaceuticalStartDateTime'))
    if decay == 'YES':
      if reference is not None:
        if abs((reference - firstFrameStart).total_seconds()) <= 1.0:
          metadata['DecayCorrection'] = 'START'
        elif injected is not None and abs((reference - injected).total_seconds()) <= 1.0:
          metadata['DecayCorrection'] = 'ADMIN'
        else:
          metadata['DecayCorrection'] = 'ENHANCED'  # Explicit non-START/ADMIN reference.
      else:
        metadata['DecayCorrection'] = 'UNKNOWN'
    elif decay == 'NO':
      metadata['DecayCorrection'] = 'NONE'

    # A single conversion factor is valid only for a common known decay reference.
    factor = None
    try:
      halfLife = float(metadata['RadionuclideHalfLife'])
      weightKg = float(metadata['PatientWeight'])
      if weightKg >= 1000:
        weightKg *= 0.001
        metadata['PatientWeight'] = str(weightKg)
      if (decay == 'YES' and reference is not None and injected is not None
          and doseBq is not None and halfLife > 0 and weightKg > 0
          and (reference - injected).total_seconds() >= -1):
        referenceDoseBq = doseBq * 2 ** (
          -(reference - injected).total_seconds() / halfLife)
        if referenceDoseBq > 0:
          candidate = weightKg * 1000.0 / referenceDoseBq
          if math.isfinite(candidate) and candidate > 0:
            factor = candidate
    except (KeyError, TypeError, ValueError, OverflowError):
      pass
    if factor is None:
      logging.warning('[dPET Enhanced] No validated SUVbw factor: verify dose, weight, '
                      'half-life, injection time and Enhanced PET decay reference')
    return metadata, factor, extra

  def examineEnhancedSingleInstance(self, files):
    """Offer an Enhanced PET loadable only for internally time-resolved files."""
    if len(files) != 1:
      return []
    sop = self._fileValue(files[0], '0008,0016')
    if sop and sop != '1.2.840.10008.5.1.4.1.1.130':
      return []
    try:
      info = self._readEnhancedFrameMetadata(files[0])
      if info is None:
        return []
      loadable = DICOMLoadable()
      loadable.files = list(files)
      desc = self._fileValue(files[0], '0008,103E') or 'Enhanced PET'
      loadable.name = f'{desc} - {info["nFrames"]} frames Enhanced PET Volume Sequence'
      loadable.tooltip = (f'Single-instance Enhanced PET, {info["nFrames"]} times x '
                          f'{info["nSlices"]} slices; units={info["units"]}. '
                          + ('Timing complete.' if all(x is not None for x in info['durations'])
                             else 'Frame durations missing: kinetic timing incomplete.'))
      loadable.selected = True
      loadable.loadAsVolumeSequence = True
      loadable.confidence = 0.99
      loadable.dPETLoadMode = 'enhanced-4d'
      return [loadable]
    except Exception as exc:
      logging.warning('[dPET Enhanced] Recognition rejected: %s', exc)
      return []

  def _iterEnhancedPixelFrames(self, filePath):
    """Stream decoded image planes; mmap native pixels on older pydicom."""
    import numpy as np
    import pydicom
    try:
      from pydicom.pixels import iter_pixels
    except ImportError:
      iter_pixels = None
    if iter_pixels is not None:
      yield from iter_pixels(filePath)
      return

    ds = pydicom.dcmread(filePath, stop_before_pixels=True)
    syntax = str(ds.file_meta.TransferSyntaxUID)
    endian = {'1.2.840.10008.1.2': '<', '1.2.840.10008.1.2.1': '<',
              '1.2.840.10008.1.2.2': '>'}.get(syntax)
    bits = int(getattr(ds, 'BitsAllocated', 0))
    if (endian is None or bits not in (8, 16, 32) or
        int(getattr(ds, 'BitsStored', 0)) != bits or
        int(getattr(ds, 'HighBit', -1)) != bits - 1):
      raise RuntimeError('Compressed/packed Enhanced PET needs pydicom 3 pixel decoding')
    deferred = pydicom.dcmread(filePath, defer_size=1024)
    raw = deferred.get_item((0x7fe0, 0x0010), keep_deferred=True)
    offset = getattr(raw, 'value_tell', None)
    length = int(ds.NumberOfFrames) * int(ds.Rows) * int(ds.Columns) * bits // 8
    if offset is None or getattr(raw, 'length', 0) < length:
      raise RuntimeError('Cannot locate complete uncompressed Enhanced PET PixelData')
    kind = 'i' if int(ds.PixelRepresentation) else 'u'
    dtype = np.dtype((endian if bits > 8 else '|') + kind + str(bits // 8))
    view = np.memmap(filePath, mode='r', offset=offset, dtype=dtype,
                     shape=(int(ds.NumberOfFrames), int(ds.Rows), int(ds.Columns)))
    for index in range(int(ds.NumberOfFrames)):
      yield view[index]

  def loadEnhancedSingleInstance(self, loadable):
    """Stream one Enhanced PET instance into a time-resolved Volume Sequence."""
    import numpy as np
    path = loadable.files[0]
    created = []
    progressbar = None
    try:
      info = self._readEnhancedFrameMetadata(path)
      if info is None:
        raise ValueError('No internal temporal dimension in Enhanced PET file')
      times, durations = info['times'], info['durations']
      if any(d is None for d in durations):
        logging.warning('[dPET Enhanced] Missing frame duration(s); kinetic timing is incomplete')
      metadata, enhancedFactor, enhancedAttributes = self._enhancedQuantitativeMetadata(
        path, times[0])
      metadata.update({'Units': info['units'], 'SUVType': info['suvType']})
      framesMeta = []
      for fi, start in enumerate(times):
        duration = durations[fi]
        framesMeta.append({
          'index': fi, 'acquisitionStartDateTime': self._formatDicomDateTime(start),
          'acquisitionStartSource': 'DICOM.FrameContent.FrameAcquisitionDateTime',
          'durationSec': duration,
          'durationSource': 'DICOM.FrameContent.FrameAcquisitionDuration' if duration is not None else 'Unavailable',
          'acquisitionEndDateTime': self._formatDicomDateTime(start + timedelta(seconds=duration)) if duration is not None else ''
        })
      kinetic = {
        'schemaVersion': 1, 'metadataSource': 'dPETImporter',
        'acquisitionKind': 'DYNAMIC', 'seriesType': 'DYNAMIC', 'wholeBody': False,
        'timingMode': 'FRAMES', 'timingComplete': all(d is not None for d in durations),
        'spatiallyVaryingTiming': False,
        'acquisitionStartDateTime': self._formatDicomDateTime(times[0]),
        'acquisitionEndDateTime': (self._formatDicomDateTime(times[-1] + timedelta(seconds=durations[-1]))
                                   if durations[-1] is not None else ''),
        'studyInstanceUID': info['studyInstanceUID'],
        'seriesInstanceUID': info['seriesInstanceUID'],
        'frameOfReferenceUID': info['frameOfReferenceUID'],
        'common': metadata, 'frames': framesMeta,
      }

      name = loadable.name
      seq = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSequenceNode',
                                                slicer.mrmlScene.GenerateUniqueName(name))
      created.append(seq)
      seq.SetIndexName('Frame')
      seq.SetIndexUnit('index')
      seq.SetAttribute('dPETImporter.LoadedBy', 'dPETImporterPlugin')
      seq.SetAttribute('dPETImporter.Source', 'DICOM')
      seq.SetAttribute('dPETImporter.Version', '0.1')
      seq.SetAttribute('dPET.FrameType', 'Enhanced-single-instance')
      seq.SetAttribute('MultiVolume.NumberOfFrames', str(info['nFrames']))
      seq.SetAttribute('MultiVolume.FrameIdentifyingDICOMTagName', 'TemporalPositionIndex')
      seq.SetAttribute('MultiVolume.FrameLabels', ','.join(str((t-times[0]).total_seconds()*1000)
                                                        for t in times))
      seq.SetAttribute('dPET.FrameTimes', json.dumps([self._formatDicomDateTime(t) for t in times]))
      seq.SetAttribute('dPET.FrameDurations', json.dumps(durations))
      # This instance owns every source plane; there is no SOP-per-slice identity.
      seq.SetAttribute('DICOM.instanceUIDs', info['sopInstanceUID'])
      for attr, field in (
          ('DICOM.StudyInstanceUID', 'studyInstanceUID'),
          ('DICOM.SeriesInstanceUID', 'seriesInstanceUID'),
          ('DICOM.FrameOfReferenceUID', 'frameOfReferenceUID'),
          ('DICOM.SOPClassUID', 'sopClassUID')):
        if info[field]:
          seq.SetAttribute(attr, info[field])
      seq.SetAttribute('dPET.Enhanced.PixelPaddingPolicy', 'NaN')
      self._applyKineticMetadata(seq, kinetic)
      for key, value in enhancedAttributes.items():
        seq.SetAttribute(key, str(value))
      suvbwFactor = enhancedFactor
      factorValid = suvbwFactor is not None and suvbwFactor > 0
      suvEnabled = settingsValue('DICOM/dPETImporterSUVEnabled', True, converter=toBool)
      doSUV = info['units'] == 'BQML' and suvEnabled and factorValid
      valueType = 'SUVbw' if info['units'] == 'GML' or doSUV else 'BQML'
      if info['units'] == 'BQML' and suvEnabled and not factorValid:
        logging.warning('[dPET Enhanced] SUVbw conversion unavailable; retaining Bq/mL')
      seq.SetAttribute('dPET.ValueType', valueType)
      seq.SetAttribute('dPET.SUVbwFactorValid', '1' if factorValid else '0')
      if factorValid:
        seq.SetAttribute('dPET.SUVbwFactor', str(suvbwFactor))

      # Source frame numbers are 1-based; spatial verification is not classic SOP-per-slice.
      provenance = {'schemaVersion': 2, 'studyInstanceUID': info['studyInstanceUID'],
                    'seriesInstanceUID': info['seriesInstanceUID'],
                    'frameOfReferenceUID': info['frameOfReferenceUID'],
                    'frames': [{'index': fi, 'instances': [{
                      'sopInstanceUID': info['sopInstanceUID'],
                      'sopClassUID': info['sopClassUID']}],
                      'sourceFrameNumbers': indices, 'spatialOrderVerified': False,
                      'spatialMapping': 'SINGLE_SOP_ENHANCED_FRAMES'}
                      for fi, indices in enumerate(info['sourceFrameNumbers'])]}
      seq.SetAttribute('dPET.DICOM.FrameReferences', json.dumps(provenance, separators=(',', ':')))
      seq.SetAttribute('dPET.DICOM.ProvenanceSchemaVersion', '2')

      frameNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLScalarVolumeNode', name + ' temporary')
      created.append(frameNode)
      transform = vtk.vtkMatrix4x4()
      for i in range(4):
        for j in range(4):
          transform.SetElement(i, j, float(info['rasMatrix'][i, j]))
      frameNode.SetIJKToRASMatrix(transform)
      progressbar = slicer.util.createProgressDialog(labelText=f'Loading {name}',
                                                    value=0, maximum=info['nFrames'],
                                                    windowModality=qt.Qt.WindowModal)
      current = None
      volume = None
      seen = 0
      paddingCount = 0
      for sourceIndex, pixels in enumerate(self._iterEnhancedPixelFrames(path)):
        if sourceIndex >= len(info['sourceToFrame']):
          raise ValueError('Decoded more image frames than listed by the DICOM header')
        temporal = info['sourceToFrame'][sourceIndex]
        if current is None or temporal != current:
          if current is not None:
            self._storeEnhancedFrame(seq, frameNode, volume, current,
                                     times, durations, valueType, suvbwFactor)
            progressbar.value = current + 1
            slicer.app.processEvents()
            if progressbar.wasCanceled:
              raise RuntimeError('Enhanced PET loading canceled')
          current = temporal
          volume = np.empty((info['nSlices'], info['rows'], info['cols']), dtype=np.float32)
        pixels = np.asarray(pixels)
        if pixels.shape != (info['rows'], info['cols']):
          raise ValueError(f'Decoded image shape differs at source frame {sourceIndex + 1}')
        units, slope, intercept, lo, hi = info['mappings'][sourceIndex]
        padding = np.zeros(pixels.shape, dtype=bool)
        if info['paddingValue'] is not None:
          first = info['paddingValue']
          last = info['paddingLimit'] if info['paddingLimit'] is not None else first
          padding = (pixels >= min(first, last)) & (pixels <= max(first, last))
        validPixels = pixels[~padding]
        if (validPixels.size and
            ((lo is not None and np.min(validPixels) < lo)
             or (hi is not None and np.max(validPixels) > hi))):
          raise ValueError(f'Pixels outside real-world mapping range at source frame {sourceIndex + 1}')
        # Large per-slice intercepts cause cancellation in float32 near zero.
        volSlice = pixels.astype(np.float64) * slope + intercept
        if doSUV:
          volSlice *= float(suvbwFactor)
        if not np.isfinite(volSlice[~padding]).all():
          raise ValueError(f'Non-finite calibrated voxels in frame {sourceIndex + 1}')
        if np.any(padding):
          paddingCount += int(np.count_nonzero(padding))
          volSlice[padding] = np.nan  # DICOM padding is not measured zero activity.
        volume[info['sourceToK'][sourceIndex]] = volSlice.astype(np.float32)
        seen += 1
        if sourceIndex % 64 == 0:
          slicer.app.processEvents()
          if progressbar.wasCanceled:
            raise RuntimeError('Enhanced PET loading canceled')
      if seen != len(info['sourceToFrame']):
        raise ValueError('Decoded fewer frames than specified by Enhanced PET metadata')
      self._storeEnhancedFrame(seq, frameNode, volume, current,
                               times, durations, valueType, suvbwFactor)
      progressbar.value = info['nFrames']
      seq.SetAttribute('dPET.Enhanced.PixelPaddingVoxelCount', str(paddingCount))
      if paddingCount:
        logging.info('[dPET Enhanced] Preserved %d DICOM padding voxels as NaN', paddingCount)

      browser = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSequenceBrowserNode', name + ' browser')
      created.append(browser)
      browser.SetAndObserveMasterSequenceNodeID(seq.GetID())
      browser.SetSaveChanges(seq, True)
      browser.SetOverwriteProxyName(seq, False)
      browser.SetAttribute('dPETImporter.LoadedBy', 'dPETImporterPlugin')
      browser.SetAttribute('dPETImporter.SequenceNodeID', seq.GetID())
      browser.SetAttribute('dPET.ValueType', valueType)
      if slicer.modules.sequences.widgetRepresentation():
        slicer.modules.sequences.widgetRepresentation().setActiveBrowserNode(browser)
      proxy = browser.GetProxyNode(seq)
      if proxy:
        created.append(proxy)
        proxy.SetName(name)
        proxy.SetAttribute('dPETImporter.LoadedBy', 'dPETImporterPlugin')
        proxy.SetAttribute('dPET.ValueType', valueType)
        for attr in ('DICOM.instanceUIDs', 'DICOM.StudyInstanceUID',
                     'DICOM.SeriesInstanceUID', 'DICOM.FrameOfReferenceUID',
                     'DICOM.SOPClassUID', 'dPET.Enhanced.PixelPaddingPolicy',
                     'dPET.Enhanced.PixelPaddingVoxelCount'):
          value = seq.GetAttribute(attr)
          if value is not None:
            proxy.SetAttribute(attr, value)
        for attr in ('dPET.FrameType', 'dPET.FrameTimes', 'dPET.FrameDurations',
                     'MultiVolume.NumberOfFrames'):
          proxy.SetAttribute(attr, seq.GetAttribute(attr))
        self._applyKineticMetadata(proxy, kinetic)
        for key, value in enhancedAttributes.items():
          proxy.SetAttribute(key, str(value))
        for attr in ('CorrectedImage', 'DecayCorrection', 'PatientWeight',
                     'RadionuclideHalfLife', 'RadionuclideTotalDose',
                     'RadionuclideStartDateTime', 'RadiopharmaceuticalStartDateTime',
                     'dPET.InjectionDateTimeSource', 'dPET.InjectionToAcquisitionOffsetSec'):
          value = seq.GetAttribute(attr)
          if value is not None:
            proxy.SetAttribute(attr, value)
        proxy.SetAttribute('dPET.SUVbwFactorValid', '1' if factorValid else '0')
        if factorValid:
          proxy.SetAttribute('dPET.SUVbwFactor', str(suvbwFactor))
        try:
          self._setProxyQuantityUnits(proxy, valueType)
          proxy.CreateDefaultDisplayNodes()
          self.setPetDicomLUT(proxy)
        except Exception:
          logging.warning('[dPET Enhanced] Display/units decoration failed', exc_info=True)
          proxy.SetAttribute('dPET.ValueType', valueType)
        selection = slicer.app.applicationLogic().GetSelectionNode()
        selection.SetReferenceActiveVolumeID(proxy.GetID())
        slicer.app.applicationLogic().PropagateVolumeSelection()
      self.addSeriesInSubjectHierarchy(loadable, proxy or seq)
      _dpetWatchPETSequenceDisplay(browser, seq, proxy)
      slicer.mrmlScene.RemoveNode(frameNode)
      created.remove(frameNode)
      return seq
    except Exception:
      logging.exception('[dPET Enhanced] Loading failed')
      for node in reversed(created):
        if node and node.GetScene():
          slicer.mrmlScene.RemoveNode(node)
      return None
    finally:
      if progressbar is not None:
        progressbar.close()

  def _storeEnhancedFrame(self, seq, frameNode, pixels, frameIndex,
                          times, durations, valueType, suvbwFactor):
    """Copy a decoded volume into the sequence with per-frame metadata."""
    slicer.util.updateVolumeFromArray(frameNode, pixels)
    if frameIndex < 4 or settingsValue('DICOM/detailedLogging', False, converter=toBool):
      import numpy as np
      valid = pixels[np.isfinite(pixels)]
      if valid.size:
        logging.info('[dPET Enhanced] Temporal frame %d voxel range: %.8g .. %.8g, '
                     'nonzero=%d/%d', frameIndex + 1, float(np.min(valid)),
                     float(np.max(valid)), int(np.count_nonzero(valid)), pixels.size)
      else:
        logging.info('[dPET Enhanced] Temporal frame %d contains only padding', frameIndex + 1)
    idx = str(frameIndex)
    timeText = self._formatDicomDateTime(times[frameIndex])
    frameNode.SetAttribute('dPETImporter.LoadedBy', 'dPETImporterPlugin')
    frameNode.SetAttribute('dPET.ValueType', valueType)
    # Preserve source identity and quantitative/correction metadata on each frame.
    sharedAttrs = ('Units', 'SUVType', 'CorrectedImage', 'DecayCorrection',
                   'PatientWeight', 'RadionuclideHalfLife', 'RadionuclideTotalDose',
                   'RadionuclideStartDateTime', 'RadiopharmaceuticalStartDateTime',
                   'RadiopharmaceuticalStartTime', 'dPET.DecayCorrectionDateTime',
                   'dPET.FirstFrameAcquisitionDateTime',
                   'dPET.InjectionDateTimeSource', 'dPET.InjectionToAcquisitionOffsetSec',
                   'DICOM.instanceUIDs', 'DICOM.StudyInstanceUID',
                   'DICOM.SeriesInstanceUID', 'DICOM.FrameOfReferenceUID',
                   'DICOM.SOPClassUID', 'dPET.Enhanced.PixelPaddingPolicy')
    sharedAttrs += tuple(attr for attr in seq.GetAttributeNames()
                         if attr.startswith('dPET.Enhanced.'))
    for attr in sharedAttrs:
      value = seq.GetAttribute(attr)
      if value is not None:
        frameNode.SetAttribute(attr, value)
    frameNode.SetAttribute('dPET.AcquisitionTime', timeText)
    frameNode.SetAttribute('dPET.SUVbwFactorValid', '1' if suvbwFactor else '0')
    if suvbwFactor:
      frameNode.SetAttribute('dPET.SUVbwFactor', str(suvbwFactor))
    else:
      frameNode.SetAttribute('dPET.SUVbwFactor', None)
    duration = durations[frameIndex]
    frameNode.SetAttribute('dPET.Duration', str(duration) if duration is not None else None)
    seq.SetAttribute(f'dPET.Frame.{idx}.AcquisitionTime', timeText)
    if duration is not None:
      seq.SetAttribute(f'dPET.Frame.{idx}.Duration', str(duration))
    seq.SetAttribute(f'dPET.Frame.{idx}.SUVbwFactorValid', '1' if suvbwFactor else '0')
    if suvbwFactor:
      seq.SetAttribute(f'dPET.Frame.{idx}.SUVbwFactor', str(suvbwFactor))
    seq.SetDataNodeAtValue(frameNode, idx)

  def examine(self,fileLists):
    dynamicEnabled = settingsValue("DICOM/dPETImporterEnabled", True, converter=toBool)
    staticEnabled = settingsValue("DICOM/dPETImporterStaticEnabled", True, converter=toBool)
    if not dynamicEnabled and not staticEnabled:
      return []

    self.detailedLogging = settingsValue('DICOM/detailedLogging', False, converter=toBool)
    dynamicLoadables = []
    enhancedLoadables = []
    staticLoadables = []
    allfiles = []
    lastFiles = []
    for files in fileLists:
      petFiles = self._petImageFiles(files)
      if not petFiles:
        continue
      enhanced = []
      if dynamicEnabled:
        enhanced = self.examineEnhancedSingleInstance(petFiles)
        if enhanced:
          enhancedLoadables += enhanced
        else:
          dynamicLoadables += self.examineFiles(petFiles)
      if staticEnabled:
        staticLoadables += self.examineStaticFiles(petFiles)
      if not enhanced:
        lastFiles = petFiles
        allfiles += petFiles

    if dynamicEnabled and (not dynamicLoadables) and lastFiles and len(allfiles) > len(lastFiles):
      dynamicLoadables += self.examineFilesMultiseries(allfiles)

    # --- annotate each dynamic loadable (MV and Sequence) with 2D/3D + SUV parsing status ---
    def annotate(loadable, isSequence=False):
      mv = getattr(loadable, "multivolume", None)
      frameType = mv.GetAttribute('dPET.FrameType') if mv else None
      if frameType == "2D-slices-per-frame":
        ft = "2D-per-frame"
      elif frameType == "3D-per-frame":
        ft = "3D-per-frame"
      else:
        ft = "unknown-dim"

      missing = self._missingSUVKeys(mv)
      suvTxt = "SUV:ok" if len(missing) == 0 else "SUV:missing"
      suffix = f" [{ft}; {suvTxt}]"
      loadable.name = (loadable.name or '') + suffix
      if loadable.tooltip:
        loadable.tooltip = loadable.tooltip + suffix
      else:
        loadable.tooltip = suffix.strip()
      if missing:
        loadable.tooltip += f" (missing: {', '.join(missing)})"

    seqLoadables = []
    if hasattr(slicer.modules, 'sequences'):
      for loadable in dynamicLoadables:
        seqL = DICOMLoadable()
        seqL.files = loadable.files
        seqL.multivolume = loadable.multivolume
        seqL.selected = loadable.selected
        seqL.confidence = loadable.confidence
        seqL.loadAsVolumeSequence = True
        seqL.dPETLoadMode = 'dynamic'
        seqL.tooltip = (loadable.tooltip or '').replace('MultiVolume', 'Volume Sequence')
        seqL.name = (loadable.name or '').replace('MultiVolume', 'Volume Sequence')
        annotate(seqL, isSequence=True)
        seqLoadables.append(seqL)

    return seqLoadables + enhancedLoadables + staticLoadables

  def nameTooltipFromFile(self, dicomFilePath, nFrames, tagName, descriptionLevel=None):
    seriesNumber = slicer.dicomDatabase.fileValue(dicomFilePath, self.tags['seriesNumber'])
    modality = slicer.dicomDatabase.fileValue(dicomFilePath, self.tags['modality'])
    if descriptionLevel=="study":
      description = slicer.dicomDatabase.fileValue(dicomFilePath, self.tags['studyDescription'])
    else:
      description = slicer.dicomDatabase.fileValue(dicomFilePath, self.tags['seriesDescription'])
    name = ''
    if seriesNumber:
      name = f'{seriesNumber}:'
    if modality:
      name = f'{name} {modality}'
    if description:
      name = f'{name} {description}'
    name = f'{name} - {nFrames} frames MultiVolume by {tagName}'
    tooltip = name
    return name, tooltip

  def examineFilesMultiseries(self,files):
    files = self._petImageFiles(files)
    if not files:
      return []
    mvNodes = self.initMultiVolumes(files,prescribedTags=['SeriesTime','AcquisitionTime'])
    loadables = []
    for mvNode in mvNodes:
      tagName = mvNode.GetAttribute('MultiVolume.FrameIdentifyingDICOMTagName')
      orderedFiles = mvNode.GetAttribute('MultiVolume.FrameFileList').split(',')
      if not self.isFrameOriginConsistent(orderedFiles, mvNode):
        continue
      loadable = DICOMLoadable()
      loadable.files = orderedFiles
      loadable.name, loadable.tooltip = self.nameTooltipFromFile(loadable.files[0], mvNode.GetNumberOfFrames(), tagName, descriptionLevel='study')
      # desc = slicer.dicomDatabase.fileValue(orderedFiles[0],self.tags['studyDescription'])
      # num = slicer.dicomDatabase.fileValue(orderedFiles[0],self.tags['seriesNumber'])
      loadable.selected = True
      loadable.multivolume = mvNode
      loadable.confidence = 1.
      loadables.append(loadable)
    return loadables

  def examineFiles(self,files):
    files = self._petImageFiles(files)
    if not files:
      return []
    subseriesLists = {}
    for f in files:
      sid = slicer.dicomDatabase.fileValue(f, self.tags['seriesInstanceUID']) or "Unknown"
      subseriesLists.setdefault(sid, []).append(f)

    loadables = []
    for sid, filelist in subseriesLists.items():
      mvNodes = self.initMultiVolumes(filelist)
      for mvNode in mvNodes:
        tagName = mvNode.GetAttribute('MultiVolume.FrameIdentifyingDICOMTagName')
        orderedFiles = mvNode.GetAttribute('MultiVolume.FrameFileList').split(',')
        if not self.isFrameOriginConsistent(orderedFiles, mvNode):
          continue
        loadable = DICOMLoadable()
        loadable.files = filelist
        loadable.name, loadable.tooltip = self.nameTooltipFromFile(filelist[0], mvNode.GetNumberOfFrames(), tagName)
        mvNode.SetName(loadable.name)
        loadable.selected = True
        loadable.multivolume = mvNode
        loadable.confidence = 1.
        loadables.append(loadable)
    return loadables

  def tm2ms(self, tm):
    if not tm or len(tm) < 2:
      return 0
    try:
      hhmmss = tm.split('.')[0]
    except:
      hhmmss = tm
    try:
      ssfrac = float('0.' + tm.split('.')[1])
    except:
      ssfrac = 0.0
    if len(hhmmss) == 6: # HHMMSS
      sec = float(hhmmss[0:2])*3600.0 + float(hhmmss[2:4])*60.0 + float(hhmmss[4:6])
    elif len(hhmmss) == 4: # HHMM
      sec = float(hhmmss[0:2])*3600.0 + float(hhmmss[2:4])*60.0
    elif len(hhmmss) == 2:
      sec = float(hhmmss[0:2])*3600.0
    else:
      sec = 0.0
    sec += ssfrac
    return sec * 1000.0

  def _multiplyVolumeByConstant(self, volumeNode, factor):
    """Multiply voxel values by factor. Output is float to avoid overflow/rounding."""
    if volumeNode is None or factor is None:
      return False
    img = volumeNode.GetImageData()
    if img is None:
      return False

    import vtk

    shiftScale = vtk.vtkImageShiftScale()
    shiftScale.SetInputData(img)
    shiftScale.SetShift(0.0)
    shiftScale.SetScale(float(factor))
    shiftScale.SetOutputScalarTypeToFloat()
    shiftScale.ClampOverflowOff()
    shiftScale.Update()

    outImg = vtk.vtkImageData()
    outImg.DeepCopy(shiftScale.GetOutput())
    volumeNode.SetAndObserveImageData(outImg)
    volumeNode.Modified()
    return True


  def _validateUnits(self, mvNode):
    units = (mvNode.GetAttribute('Units') or '').upper()
    suvType = (mvNode.GetAttribute('SUVType') or '').upper()

    if not units:
      logging.error("[dPET] Missing Units (0054,1001)")
      return False, None

    if units == "BQML":
      return True, "BQML"

    # Not BQML → must be SUV-like
    if units == "GML":
      # SlicerDynamicPET currently supports body-weight normalized SUV only.
      if not suvType:
        logging.warning("[dPET] Units=GML but SUVType missing; assuming BW")
        suvType = "BW"
      if suvType not in ("BW", "SUVBW"):
        logging.error(f"[dPET] Unsupported SUV type: {suvType}. Only SUVbw is supported.")
        return False, None
      return True, "SUV"

    # Any other unit → unsupported
    logging.error(f"[dPET] Unsupported Units: {units}. Please follow the DICOM standard.")
    return False, None

  def initMultiVolumes(self, files, prescribedTags=None):
    files = self._petImageFiles(files)
    if not files:
      return []
    tag2ValueFileList = {}
    multivolumes = []
    consideredTags = list(self.multiVolumeTags.keys()) if prescribedTags is None else list(prescribedTags)
    # iterate over all files
    tagsToIgnore = []
    for file in files:
      for frameTag in tagsToIgnore:
        if frameTag in consideredTags:
          consideredTags.remove(frameTag)
      tagsToIgnore = []

      for frameTag in list(consideredTags):
        tagMap = tag2ValueFileList.setdefault(frameTag, {})
        tagValueStr = slicer.dicomDatabase.fileValue(file, self.tags[frameTag])
        if not tagValueStr:
          tagsToIgnore.append(frameTag)
          continue
        if frameTag in ('AcquisitionTime','SeriesTime','ContentTime'):
          tagValue = self.tm2ms(tagValueStr)
        elif frameTag == 'TemporalPositionIdentifier':
          try:
            tagValue = float(tagValueStr)
          except:
            continue
        else:
          try:
            tagValue = float(tagValueStr)
          except:
            continue
        tagMap.setdefault(tagValue, []).append(file)

    # iterate over the parsed items and decide which ones can qualify as mv
    for frameTag in consideredTags:
      tagValue2FileList = tag2ValueFileList.get(frameTag)
      if not tagValue2FileList or len(tagValue2FileList) < 2:
        continue
      tagValues = sorted(tagValue2FileList.keys())

      slicesPerFrame = {}
      for tv in tagValues:
        nSlices = len(tagValue2FileList[tv])
        slicesPerFrame.setdefault(nSlices, []).append(tv)
      if len(slicesPerFrame) > 1:
        continue

      geometryValid = True
      for tv in tagValues:
        positions = set()
        orientations = set()
        frameFileList = tagValue2FileList[tv]
        for f in frameFileList:
          positions.add(slicer.dicomDatabase.fileValue(f, self.tags['position']))
          orientations.add(slicer.dicomDatabase.fileValue(f, self.tags['orientation']))
        if len(positions) != len(frameFileList) or len(orientations) != 1:
          geometryValid = False
          break
      if not geometryValid:
        continue

      # build mvNode describing this grouping
      frameFileListStr = ",".join([f for tv in tagValues for f in tagValue2FileList[tv]])
      frameLabelsArray = vtk.vtkDoubleArray()
      base = tagValues[0]
      for tv in tagValues:
        if frameTag in ('SeriesTime','AcquisitionTime','ContentTime'):
          frameLabelsArray.InsertNextValue(tv - base)
        else:
          frameLabelsArray.InsertNextValue(tv)

      # compute per-frame times and durations (best-effort)
      frame_times = []
      frame_durations = []
      # common duration tags (vendor-dependent)
      dur_tags = ["0067,1004", "0018,1242"]  # (0067,1004) expected in seconds, (0018,1242) expected in ms
      time_tags = ["0008,0032", "0008,002A", "0008,0031", "0008,0033", "0018,1075"]  # AcquisitionTime, AcquisitionDateTime, SeriesTime, ContentTime, Private
      for tv in tagValues:
        frame_files = tagValue2FileList[tv]
        # try duration
        duration = None
        for dt in dur_tags:
          val = slicer.dicomDatabase.fileValue(frame_files[0], dt)
          if val not in (None, ""):
            # try parse numeric, allow backslash-separated values
            s = str(val)
            try:
              duration = float(s.split('\\')[0])
              duration = duration/1000 if dt=="0018,1242" else duration
              # heuristics: many vendors use ms or sec; keep raw value; caller must interpret units
            except:
              duration = None
            break
        # try acquisition datetime/time
        acqtime = None
        for tt in time_tags:
          val = slicer.dicomDatabase.fileValue(frame_files[0], tt)
          if val not in (None, ""):
            acqtime = val
            break
        frame_times.append(acqtime)
        frame_durations.append(duration)

      mvNode = slicer.mrmlScene.CreateNodeByClass('vtkMRMLMultiVolumeNode')
      mvNode.UnRegister(None)
      mvNode.SetAttribute("MultiVolume.FrameLabels", ",".join(str(x) for x in [frameLabelsArray.GetValue(i) for i in range(frameLabelsArray.GetNumberOfTuples())]))
      mvNode.SetAttribute("MultiVolume.FrameIdentifyingDICOMTagName", frameTag)
      mvNode.SetAttribute('MultiVolume.NumberOfFrames', str(len(tagValues)))
      mvNode.SetAttribute('MultiVolume.FrameIdentifyingDICOMTagUnits', self.multiVolumeTagsUnits.get(frameTag, ""))
      mvNode.SetAttribute('MultiVolume.FrameFileList', frameFileListStr)
      mvNode.SetNumberOfFrames(len(tagValues))
      mvNode.SetLabelName(self.multiVolumeTagsUnits.get(frameTag, ""))
      mvNode.SetLabelArray(frameLabelsArray)
      # store computed arrays (JSON)
      try:
        mvNode.SetAttribute('dPET.FrameTimes', json.dumps(frame_times))
        mvNode.SetAttribute('dPET.FrameDurations', json.dumps(frame_durations))
      except Exception:
        # json may fail for weird objects; fallback to string join
        mvNode.SetAttribute('dPET.FrameTimes', ",".join([str(x) for x in frame_times]))
        mvNode.SetAttribute('dPET.FrameDurations', ",".join([str(x) for x in frame_durations]))
      # determine type: 3D-per-frame if 1 file per frame, else 2D-slices-per-frame
      nFilesTotal = sum(len(tagValue2FileList[tv]) for tv in tagValues)
      files_per_frame = int(nFilesTotal / len(tagValues)) if len(tagValues)>0 else 0
      frame_type = "3D-per-frame" if files_per_frame == 1 else "2D-slices-per-frame"
      mvNode.SetAttribute('dPET.FrameType', frame_type)

      # Determine first-frame acquisition datetime explicitly
      firstFrameFile = tagValue2FileList[tagValues[0]][0]

      acqDateTime = slicer.dicomDatabase.fileValue(firstFrameFile, "0008,002A")  # AcquisitionDateTime
      firstFrameDt = None

      if acqDateTime not in (None, ""):
        firstFrameDt = _parse_dicom_datetime_value(acqDateTime)
      else:
        firstFrameDate = (
          slicer.dicomDatabase.fileValue(firstFrameFile, "0008,0022") or  # AcquisitionDate
          slicer.dicomDatabase.fileValue(firstFrameFile, "0008,0021") or  # SeriesDate
          slicer.dicomDatabase.fileValue(firstFrameFile, "0008,0020")     # StudyDate
        )
        firstFrameTime = (
          slicer.dicomDatabase.fileValue(firstFrameFile, "0008,0032") or  # AcquisitionTime
          slicer.dicomDatabase.fileValue(firstFrameFile, "0008,0031") or  # SeriesTime
          slicer.dicomDatabase.fileValue(firstFrameFile, "0008,0033")     # ContentTime
        )
        firstFrameDt = _parse_dicom_datetime(firstFrameDate, firstFrameTime)

      if firstFrameDt is not None:
        mvNode.SetAttribute('dPET.FirstFrameAcquisitionDateTime', firstFrameDt.strftime("%Y%m%d%H%M%S"))

      def setAttrIfFound(attrName, dicomTag):
        v = slicer.dicomDatabase.fileValue(firstFrameFile, dicomTag)
        if v not in (None, ""):
          mvNode.SetAttribute(attrName, v)
        else:
          tagres = self.find_all_dicom_tag_values(firstFrameFile, dicomTag, return_element=False)
          if len(tagres)>0:
            v = tagres[0]
            mvNode.SetAttribute(attrName, str(v))

      # Patient
      setAttrIfFound('PatientWeight', "0010,1030")  # Patient's Weight (kg)
      setAttrIfFound('Units', "0054,1001")     # Units
      setAttrIfFound('SUVType', "0054,1006")   # SUV Type (optional)

      # PET corrections
      setAttrIfFound('DecayCorrection', "0054,1102")   # Decay Correction
      setAttrIfFound('CorrectedImage', "0028,0051")    # Corrected Image (often contains ATTN/DECY/DECAY)

      setAttrIfFound('RadionuclideHalfLife', "0018,1075")      # Radionuclide Half Life
      setAttrIfFound('RadionuclideTotalDose', "0018,1074")     # Radionuclide Total Dose
      setAttrIfFound('RadiopharmaceuticalStartDateTime', "0018,1078")  # DICOM DT, preferred
      setAttrIfFound('RadiopharmaceuticalStartTime', "0018,1072")      # DICOM TM fallback

      needNested = (
        not mvNode.GetAttribute('RadionuclideHalfLife') or
        not mvNode.GetAttribute('RadionuclideTotalDose') or
        (not mvNode.GetAttribute('RadiopharmaceuticalStartDateTime') and
         not mvNode.GetAttribute('RadiopharmaceuticalStartTime'))
      )
      if needNested:
        nested = self._getRadiopharmNested(firstFrameFile)
        for k, v in nested.items():
          if v and not mvNode.GetAttribute(k):
            mvNode.SetAttribute(k, v)

      setAttrIfFound('StudyDate', "0008,0020")
      setAttrIfFound('SeriesDate', "0008,0021")
      setAttrIfFound('SeriesTime', "0008,0031")

      resolvedStartDt, startSource = self._resolveRadiopharmaceuticalStartDateTime(
        mvNode, firstFrameDt)
      if resolvedStartDt is not None:
        resolvedText = resolvedStartDt.strftime("%Y%m%d%H%M%S")
        # Keep the exact DICOM attribute separate when it exists, while
        # retaining this backwards-compatible normalized value for consumers.
        mvNode.SetAttribute('RadionuclideStartDateTime', resolvedText)
        mvNode.SetAttribute('dPET.InjectionDateTimeSource', startSource)

        if firstFrameDt is not None:
          rawOffsetSec = (firstFrameDt - resolvedStartDt).total_seconds()
          mvNode.SetAttribute(
            'dPET.InjectionToAcquisitionOffsetSec',
            str(float(rawOffsetSec)))
      else:
        mvNode.SetAttribute('dPET.InjectionDateTimeSource', 'Unavailable')

      multivolumes.append(mvNode)

    return multivolumes

  def find_all_dicom_tag_values(self, source, tag, return_element=False):
      import pydicom
      # Normalize tag
      if isinstance(tag, str):
          group, element = tag.split(',')
          tag = pydicom.tag.Tag(int(group, 16), int(element, 16))
      else:
          tag = pydicom.tag.Tag(tag)

      # Load dataset if needed
      if isinstance(source, str):
          ds = pydicom.dcmread(source, stop_before_pixels=True, force=True)
      else:
          ds = source

      results = []

      def recursive_search(dataset):
          for elem in dataset:
              if elem.tag == tag:
                  results.append(elem if return_element else elem.value)

              if elem.VR == "SQ":  # Sequence
                  for item in elem.value:
                      recursive_search(item)

      recursive_search(ds)
      return results

  def isFrameOriginConsistent(self, files, mvNode):
    nFrames = mvNode.GetNumberOfFrames()
    nFiles = len(files)
    if nFrames == 0 or nFiles % nFrames != 0:
      return False
    filesPerFrame = int(nFiles / nFrames)
    scalarVolumePlugin = slicer.modules.dicomPlugins['DICOMScalarVolumePlugin']()
    firstOrigin = None
    for f in range(nFrames):
      frameFiles = files[f*filesPerFrame:(f+1)*filesPerFrame]
      svs = scalarVolumePlugin.examine([frameFiles])
      if not svs:
        return False
      pos = None
      if len(svs[0].files)>1:
        pos = slicer.dicomDatabase.fileValue(svs[0].files[0], self.tags['position'])
      else:
        all_pos = self.find_all_dicom_tag_values(svs[0].files[0], self.tags['position'], return_element=False)
        pos = all_pos[0]

      if not pos:
        return False
      if len(pos)==3:
        origin = [float(x) for x in pos]
      else:
        origin = [float(x) for x in pos.split('\\')]
      if firstOrigin is None:
        firstOrigin = origin
      else:
        for a,b in zip(firstOrigin, origin):
          if abs(a-b) > self.epsilon:
            return False
    return True

  def setPetDicomLUT(self, volumeNode):
    """Apply the standard Slicer PET display to a loaded PET volume.

    The display is made deterministic for both static PET and dynamic PET
    proxies:
      * use the DICOM-standard PET procedural color scale;
      * initialize a valid W/L from the final voxel scalar range;
      * leave Auto Window/Level ON;
      * enable interpolation.

    Initializing W/L before turning Auto back on avoids a stale display state
    inherited from the generic scalar-volume loader while preserving Slicer's
    normal automatic W/L behavior afterwards.
    """
    if not volumeNode:
      return False

    displayNode = volumeNode.GetDisplayNode()
    if displayNode is None:
      volumeNode.CreateDefaultDisplayNodes()
      displayNode = volumeNode.GetDisplayNode()
    if displayNode is None:
      logging.error("[dPET] Could not create PET volume display node.")
      return False

    # Prefer Slicer's built-in PET-DICOM procedural color node.
    petColorNode = slicer.mrmlScene.GetFirstNodeByName("PET-DICOM")

    # Be defensive across Slicer builds where the singleton node may not yet
    # have been instantiated in the scene.
    if petColorNode is None:
      try:
        colorLogicClass = getattr(slicer, "vtkMRMLColorLogic", None)
        petClass = getattr(slicer, "vtkMRMLPETProceduralColorNode", None)
        if colorLogicClass is not None and petClass is not None:
          petColorNodeID = colorLogicClass.GetPETColorNodeID(petClass.PETDICOM)
          if petColorNodeID:
            petColorNode = slicer.mrmlScene.GetNodeByID(petColorNodeID)
      except Exception:
        petColorNode = None

    # Final fallback: create a PET procedural color node explicitly and set it
    # to the DICOM-standard PET palette.
    if petColorNode is None:
      try:
        petColorNode = slicer.mrmlScene.AddNewNodeByClass(
          "vtkMRMLPETProceduralColorNode", "PET-DICOM")
        if hasattr(petColorNode, "SetTypeToDICOM"):
          petColorNode.SetTypeToDICOM()
        elif hasattr(petColorNode, "SetType"):
          petClass = getattr(slicer, "vtkMRMLPETProceduralColorNode", None)
          if petClass is not None:
            petColorNode.SetType(petClass.PETDICOM)
      except Exception as error:
        logging.error(f"[dPET] Could not create PET-DICOM color scale: {error}")
        petColorNode = None

    if petColorNode is not None:
      displayNode.SetAndObserveColorNodeID(petColorNode.GetID())
    else:
      logging.error("[dPET] PET-DICOM color node could not be resolved.")

    imageData = volumeNode.GetImageData()
    if imageData is not None:
      try:
        scalarMin, scalarMax = imageData.GetScalarRange()
        scalarMin = float(scalarMin)
        scalarMax = float(scalarMax)
        if scalarMax > scalarMin:
          # Seed a valid display range first. SetWindowLevelMinMax does not
          # represent the desired final state; Auto W/L is explicitly turned
          # back on immediately afterwards.
          displayNode.AutoWindowLevelOff()
          displayNode.SetWindowLevelMinMax(scalarMin, scalarMax)
      except Exception as error:
        logging.warning(f"[dPET] PET scalar-range W/L initialization failed: {error}")

    # Leave the node explicitly in automatic W/L mode. Use the convenience
    # On() methods because they also make the intended MRML state unambiguous.
    displayNode.AutoWindowLevelOn()
    displayNode.InterpolateOn()

    # Trigger the display pipeline after the final image values and color node
    # are both in place.
    if imageData is not None:
      imageData.Modified()
    displayNode.Modified()
    volumeNode.Modified()

    # Store simple diagnostics to make runtime verification easy.
    volumeNode.SetAttribute(
      "dPET.Display.AutoWindowLevel",
      "1" if displayNode.GetAutoWindowLevel() else "0")
    volumeNode.SetAttribute(
      "dPET.Display.ColorNode",
      str(displayNode.GetColorNodeID() or ""))

    return bool(
      displayNode.GetAutoWindowLevel()
      and petColorNode is not None
      and displayNode.GetColorNodeID() == petColorNode.GetID())


  def _setProxyQuantityUnits(self, volumeNode, valueType):
    if volumeNode is None:
      return

    # These MRML APIs expect vtkCodedEntry* (VTK object), not Python strings.
    if not (hasattr(volumeNode, "SetVoxelValueQuantity") and hasattr(volumeNode, "SetVoxelValueUnits")):
      volumeNode.SetAttribute("dPET.ValueType", valueType)
      return

    # vtkCodedEntry is in Slicer MRML; depending on Slicer/VTK build it may be exposed via vtk or slicer.
    def makeCodedEntry(codeValue, codingScheme, codeMeaning):
      CodedEntryClass = getattr(vtk, "vtkCodedEntry", None) or getattr(slicer, "vtkCodedEntry", None)
      if CodedEntryClass is None:
        raise RuntimeError("vtkCodedEntry class is not available in this Slicer build.")
      ce = CodedEntryClass()
      # Convenience API supported by vtkCodedEntry
      ce.SetValueSchemeMeaning(str(codeValue), str(codingScheme), str(codeMeaning))
      return ce

    if valueType == "SUVbw":
      # Quantity: no universal “standard” code is guaranteed here; use a private scheme
      quantity = makeCodedEntry("SUVbw", "99dPET", "Standardized Uptake Value (body weight)")

      # Units: SUV is dimensionless -> UCUM “1”
      units = makeCodedEntry("1", "UCUM", "no units")
    else:
      # If you're keeping the image in BQML, your earlier label says "Bq/mL"
      quantity = makeCodedEntry("Radioactivity concentration", "99dPET", "Radioactivity concentration")
      units = makeCodedEntry("Bq/mL", "UCUM", "Becquerel per milliliter")

    volumeNode.SetVoxelValueQuantity(quantity)
    volumeNode.SetVoxelValueUnits(units)
    volumeNode.SetAttribute("dPET.ValueType", valueType)

  @staticmethod
  def _parseDICOMVector(value, expectedLength=None):
    """Parse a DICOM multi-value string into floats without re-reading files."""
    if value in (None, ""):
      return None
    try:
      if isinstance(value, (list, tuple)):
        values = [float(item) for item in value]
      else:
        values = [float(item) for item in str(value).split('\\') if item != '']
      if expectedLength is not None and len(values) != expectedLength:
        return None
      return values
    except Exception:
      return None

  def _canonicalizeLoadedFrameProvenance(self, frameNode, frameRecord, positionToleranceMm=1e-2):
    """Put classic single-frame source instances in loaded Slicer K order.

    DICOMScalarVolumePlugin sets DICOM.instanceUIDs so that UID i corresponds
    to loaded voxel slice K=i, including any slice-order reversal performed to
    obtain a right-handed IJK coordinate system.  Use that as the canonical
    ordering, then independently verify the mapping against the Image Position
    (Patient) values captured from the indexed source metadata.

    Returns True when source identity is consistent.  For classic slice stacks
    ``spatialOrderVerified`` is True only after the geometric position check.
    A single-SOP 3D/multiframe source is retained as valid source identity but
    is explicitly marked as not having a per-K SOP mapping.
    """
    if frameNode is None or frameRecord is None:
      return False

    uidText = frameNode.GetAttribute('DICOM.instanceUIDs') or ''
    orderedUIDs = [uid for uid in uidText.split() if uid]
    instances = list(frameRecord.get('instances') or [])
    if not orderedUIDs or not instances:
      return False

    byUID = {}
    for instance in instances:
      uid = str(instance.get('sopInstanceUID') or '')
      if not uid or uid in byUID:
        return False
      byUID[uid] = instance

    if len(orderedUIDs) != len(instances) or set(orderedUIDs) != set(byUID):
      logging.warning(
        '[dPET] Loaded DICOM.instanceUIDs do not match the persisted source '
        'instance set; exact spatial provenance was not accepted.')
      frameRecord['spatialOrderVerified'] = False
      return False

    orderedInstances = [byUID[uid] for uid in orderedUIDs]
    frameRecord['instances'] = orderedInstances
    frameRecord['loadedInstanceUIDOrder'] = orderedUIDs

    imageData = frameNode.GetImageData()
    if imageData is None:
      frameRecord['spatialOrderVerified'] = False
      return False

    numberOfSlices = int(imageData.GetDimensions()[2])

    # One SOP owning the complete 3D/multiframe volume is a valid source
    # identity, but there is no one-SOP-per-K mapping to verify here.
    if len(orderedUIDs) == 1 and numberOfSlices > 1:
      frameRecord['spatialMapping'] = 'SINGLE_SOP_VOLUME'
      frameRecord['spatialOrderVerified'] = False
      return True

    if len(orderedUIDs) != numberOfSlices:
      logging.warning(
        '[dPET] Source instance count does not match loaded PET K dimension; '
        'exact per-slice DICOM provenance was not accepted.')
      frameRecord['spatialOrderVerified'] = False
      return False

    ijkToRAS = vtk.vtkMatrix4x4()
    frameNode.GetIJKToRASMatrix(ijkToRAS)

    # Verify UID[k] against the physical position of loaded Slicer slice K=k.
    for k, instance in enumerate(orderedInstances):
      position = instance.get('imagePositionPatient')
      if not isinstance(position, (list, tuple)) or len(position) != 3:
        frameRecord['spatialOrderVerified'] = False
        return False

      ras = [0.0, 0.0, 0.0, 1.0]
      ijkToRAS.MultiplyPoint([0.0, 0.0, float(k), 1.0], ras)
      loadedLPS = [-float(ras[0]), -float(ras[1]), float(ras[2])]

      if max(abs(float(position[i]) - loadedLPS[i]) for i in range(3)) > float(positionToleranceMm):
        logging.warning(
          f'[dPET] Source UID {orderedUIDs[k]} position does not match loaded '
          f'Slicer slice K={k}; exact spatial provenance was not accepted.')
        frameRecord['spatialOrderVerified'] = False
        return False

    # Verify that the source in-plane orientation agrees with the loaded
    # Slicer I/J axes after RAS -> DICOM LPS conversion.
    def normalizedLPSAxis(column):
      vector = [
        -float(ijkToRAS.GetElement(0, column)),
        -float(ijkToRAS.GetElement(1, column)),
         float(ijkToRAS.GetElement(2, column)),
      ]
      norm = sum(component * component for component in vector) ** 0.5
      return [component / norm for component in vector] if norm > 0.0 else None

    iAxisLPS = normalizedLPSAxis(0)
    jAxisLPS = normalizedLPSAxis(1)
    if iAxisLPS is None or jAxisLPS is None:
      frameRecord['spatialOrderVerified'] = False
      return False

    expectedOrientation = iAxisLPS + jAxisLPS
    for instance in orderedInstances:
      orientation = instance.get('imageOrientationPatient')
      if not isinstance(orientation, (list, tuple)) or len(orientation) != 6:
        frameRecord['spatialOrderVerified'] = False
        return False
      if max(abs(float(orientation[i]) - expectedOrientation[i]) for i in range(6)) > 1e-5:
        logging.warning(
          '[dPET] Source ImageOrientationPatient does not match the loaded PET '
          'I/J axes; exact spatial provenance was not accepted.')
        frameRecord['spatialOrderVerified'] = False
        return False

    frameRecord['spatialMapping'] = 'CLASSIC_SOP_PER_K'
    frameRecord['spatialOrderVerified'] = True
    return True


  def _buildDynamicRTExportProvenance(self, files, nFrames, filesPerFrame):
    """Build compact frame-wise DICOM provenance from the indexed DICOM DB.

    No source file is opened here. The DICOM database has already indexed the
    tags used by dPETImporter, so retaining them while loading adds negligible
    work compared with creating the PET volumes themselves.

    One temporal frame owns one or more source SOP instances:
      * 3D-per-frame: one instance;
      * 2D-slices-per-frame: one instance per spatial slice.
    """
    if not files or nFrames < 1 or filesPerFrame < 1:
      return None

    db = slicer.dicomDatabase
    firstFile = files[0]

    def value(filePath, tag):
      result = db.fileValue(filePath, tag)
      return str(result) if result not in (None, "") else ""

    patientStudyTags = {
      'SpecificCharacterSet': '0008,0005',
      'PatientName': '0010,0010',
      'PatientID': '0010,0020',
      'PatientBirthDate': '0010,0030',
      'PatientSex': '0010,0040',
      'PatientAge': '0010,1010',
      'PatientSize': '0010,1020',
      'PatientWeight': '0010,1030',
      'StudyDate': '0008,0020',
      'StudyTime': '0008,0030',
      'AccessionNumber': '0008,0050',
      'StudyID': '0020,0010',
      'ReferringPhysicianName': '0008,0090',
      'PerformingPhysicianName': '0008,1050',
      'OperatorsName': '0008,1070',
      'InstitutionName': '0008,0080',
      'InstitutionAddress': '0008,0081',
    }

    provenance = {
      'schemaVersion': 2,
      'studyInstanceUID': value(firstFile, '0020,000D'),
      'seriesInstanceUID': value(firstFile, '0020,000E'),
      'frameOfReferenceUID': value(firstFile, '0020,0052'),
      'patientStudy': {},
      'frames': [],
    }
    for keyword, tag in patientStudyTags.items():
      tagValue = value(firstFile, tag)
      if tagValue:
        provenance['patientStudy'][keyword] = tagValue

    for frameIndex in range(nFrames):
      frameFiles = files[frameIndex * filesPerFrame:(frameIndex + 1) * filesPerFrame]
      instances = []
      for filePath in frameFiles:
        sopInstanceUID = value(filePath, '0008,0018')
        sopClassUID = value(filePath, '0008,0016')
        if not sopInstanceUID or not sopClassUID:
          # A frame without stable DICOM identity is not useful provenance.
          return None
        instance = {
          'sopInstanceUID': sopInstanceUID,
          'sopClassUID': sopClassUID,
        }
        position = self._parseDICOMVector(value(filePath, '0020,0032'), 3)
        orientation = self._parseDICOMVector(value(filePath, '0020,0037'), 6)
        if position is not None:
          instance['imagePositionPatient'] = position
        if orientation is not None:
          instance['imageOrientationPatient'] = orientation
        instanceNumber = value(filePath, '0020,0013')
        if instanceNumber:
          instance['instanceNumber'] = instanceNumber

        instances.append(instance)

      provenance['frames'].append({
        'index': frameIndex,
        'instances': instances,
      })

    if (not provenance['studyInstanceUID']
        or not provenance['seriesInstanceUID']
        or not provenance['frameOfReferenceUID']):
      return None
    return provenance


  def _loadStaticPET(self, loadable):
    """Load a STATIC/WHOLE BODY PET as scalar volume and persist kinetic metadata."""
    scalarPluginClass = slicer.modules.dicomPlugins.get('DICOMScalarVolumePlugin')
    if scalarPluginClass is None:
      logging.error('[dPET static] DICOMScalarVolumePlugin is unavailable.')
      return None
    scalarPlugin = scalarPluginClass()

    scalarLoadable = getattr(loadable, 'dPETScalarLoadable', None)
    if scalarLoadable is None:
      candidates = scalarPlugin.examine([loadable.files]) or []
      candidates = [item for item in candidates if getattr(item, 'files', None)]
      if not candidates:
        logging.error('[dPET static] Scalar PET interpretation failed.')
        return None
      scalarLoadable = max(candidates, key=lambda item: len(item.files))
    scalarLoadable.name = loadable.name.replace(' [dPET kinetic metadata]', '')

    volumeNode = scalarPlugin.load(scalarLoadable)
    if volumeNode is None:
      logging.error('[dPET static] Scalar PET load failed.')
      return None

    metadata = getattr(loadable, 'dPETKineticMetadata', None)
    if not metadata:
      metadata = self._buildStaticKineticMetadata(scalarLoadable.files)
    self._applyKineticMetadata(volumeNode, metadata)
    volumeNode.SetAttribute('dPETImporter.LoadedBy', 'dPETImporterPlugin')
    volumeNode.SetAttribute('dPETImporter.Version', '0.2')
    volumeNode.SetAttribute('dPETImporter.Source', 'DICOM')
    volumeNode.SetAttribute('dPETImporter.StaticKineticPET', '1')

    ok, unitType = self._validateUnits(volumeNode)
    doSUV = settingsValue('DICOM/dPETImporterSUVEnabled', True, converter=toBool)
    factor = compute_suvbw_factor(volumeNode) if ok else None
    factorValid = factor is not None and factor > 0.0
    spatialTiming = bool(metadata and metadata.get('spatiallyVaryingTiming'))

    valueType = (volumeNode.GetAttribute('Units') or 'Unknown').upper()
    if ok and unitType == 'SUV':
      valueType = 'SUVbw'
    elif ok and unitType == 'BQML':
      valueType = 'BQML'
      # A single whole-volume factor is intentionally used only when every
      # source image shares the same acquisition interval. For spatially
      # varying whole-body timing, preserve Bq/mL and the metadata instead of
      # applying a potentially incorrect global decay/SUV normalization.
      if doSUV and factorValid and not spatialTiming:
        if self._multiplyVolumeByConstant(volumeNode, factor):
          valueType = 'SUVbw'
      elif doSUV and spatialTiming:
        logging.warning(
          '[dPET static] Spatially varying whole-body timing detected; keeping Bq/mL '
          'instead of applying one global SUVbw factor. Multi-timepoint analysis can '
          'use the persisted timing/quantitative metadata to validate compatibility.')

    volumeNode.SetAttribute('dPET.ValueType', valueType)
    globalFactorValid = bool(factorValid and not spatialTiming)
    if globalFactorValid:
      volumeNode.SetAttribute('dPET.SUVbwFactor', str(float(factor)))
      volumeNode.SetAttribute('dPET.SUVbwFactorValid', '1')
    else:
      # A factor computed from the earliest slice is not a validated whole-volume
      # factor when spatial acquisition timing varies. Do not advertise it as such.
      volumeNode.SetAttribute('dPET.SUVbwFactorValid', '0')

    if valueType in ('SUVbw', 'BQML'):
      self._setProxyQuantityUnits(volumeNode, valueType)
    self.setPetDicomLUT(volumeNode)
    _dpetApplyCollapsedAutoWindowFallback(volumeNode)

    appLogic = slicer.app.applicationLogic()
    if appLogic:
      selectionNode = appLogic.GetSelectionNode()
      selectionNode.SetReferenceActiveVolumeID(volumeNode.GetID())
      appLogic.PropagateVolumeSelection()

    timingMode = metadata.get('timingMode', 'INCOMPLETE') if metadata else 'INCOMPLETE'
    logging.info(
      f"[dPET static] Loaded kinetic-ready PET '{volumeNode.GetName()}' "
      f"(timing={timingMode}, valueType={valueType}).")
    return volumeNode

  def load(self,loadable):
    """Load dynamic PET as sequence or optional static PET as scalar volume."""
    candidateFiles = list(getattr(loadable, 'files', []) or [])
    if not candidateFiles or any(not self._isPETImageFile(path) for path in candidateFiles):
      return None

    if getattr(loadable, 'dPETLoadMode', None) == 'static':
      return self._loadStaticPET(loadable)
    if getattr(loadable, 'dPETLoadMode', None) == 'enhanced-4d':
      return self.loadEnhancedSingleInstance(loadable)

    try:
      mvNode = loadable.multivolume
    except AttributeError:
      return None

    ok, unitType = self._validateUnits(mvNode)
    if not ok:
      logging.error("[dPET] Invalid or unsupported PET units. Skipping load.")
      return None

    nFrames = int(mvNode.GetAttribute('MultiVolume.NumberOfFrames'))
    files = mvNode.GetAttribute('MultiVolume.FrameFileList').split(',')
    nFiles = len(files)
    if nFrames == 0 or nFiles % nFrames != 0:
      return None
    filesPerFrame = int(nFiles/nFrames)
    baseName = loadable.name

    volumeSequenceNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSequenceNode", slicer.mrmlScene.GenerateUniqueName(baseName))
    volumeSequenceNode.SetIndexName("")
    volumeSequenceNode.SetIndexUnit("")
    volumeSequenceNode.SetAttribute("dPETImporter.LoadedBy", "dPETImporterPlugin")
    volumeSequenceNode.SetAttribute("dPETImporter.Version", "0.1")  # optional
    volumeSequenceNode.SetAttribute("dPETImporter.Source", "DICOM") # optional
    for attr in mvNode.GetAttributeNames():
      volumeSequenceNode.SetAttribute(attr, mvNode.GetAttribute(attr))

    kineticMetadata = self._buildDynamicKineticMetadata(
      files, nFrames, filesPerFrame, mvNode)
    self._applyKineticMetadata(volumeSequenceNode, kineticMetadata)

    # Capture source identity/geometry metadata while the indexed DICOM series
    # is still available.  Per-frame instance ordering is finalized only after
    # DICOMScalarVolumePlugin has loaded each frame, because that plugin is the
    # authority for the final Slicer K ordering.
    dynamicRTProvenance = self._buildDynamicRTExportProvenance(
      files, nFrames, filesPerFrame)
    if dynamicRTProvenance is not None:
      volumeSequenceNode.SetAttribute(
        'dPET.DICOM.StudyInstanceUID', dynamicRTProvenance['studyInstanceUID'])
      volumeSequenceNode.SetAttribute(
        'dPET.DICOM.SeriesInstanceUID', dynamicRTProvenance['seriesInstanceUID'])
      frameOfReferenceUID = dynamicRTProvenance.get('frameOfReferenceUID') or ''
      if frameOfReferenceUID:
        volumeSequenceNode.SetAttribute(
          'dPET.DICOM.FrameOfReferenceUID', frameOfReferenceUID)
    else:
      logging.warning(
        '[dPET] Could not capture complete DICOM source provenance at import time.')

    try:
      frame_times = json.loads(mvNode.GetAttribute('dPET.FrameTimes') or '[]')
    except Exception:
      times_attr = mvNode.GetAttribute('dPET.FrameTimes') or ''
      frame_times = times_attr.split(',') if times_attr else []

    try:
      frame_durations = json.loads(mvNode.GetAttribute('dPET.FrameDurations') or '[]')
    except Exception:
      dur_attr = mvNode.GetAttribute('dPET.FrameDurations') or ''
      frame_durations = [float(x) if x else None for x in dur_attr.split(',')] if dur_attr else []

    try:
      if frame_times:
        volumeSequenceNode.SetAttribute('dPET.FrameTimes', json.dumps(frame_times))
      if frame_durations:
        volumeSequenceNode.SetAttribute('dPET.FrameDurations', json.dumps(frame_durations))
    except Exception:
      pass

    scalarVolumePlugin = slicer.modules.dicomPlugins['DICOMScalarVolumePlugin']()
    progressbar = slicer.util.createProgressDialog(
        labelText=f"Loading {baseName}",
        value=0,
        maximum=nFrames,
        windowModality=qt.Qt.WindowModal
    )
    sequenceValueType = "BQML" if unitType == "BQML" else "SUVbw"
    try:
      doSUV = settingsValue("DICOM/dPETImporterSUVEnabled", True, converter=toBool)

      # Compute the physical Bq/mL -> SUVbw factor independently of whether
      # voxel values are converted on load. This allows DynamicPET to move
      # safely between SUVbw and activity concentration later without changing
      # the stored image values.
      suvbwFactor = compute_suvbw_factor(mvNode)
      factorValid = suvbwFactor is not None and suvbwFactor > 0

      sequenceSUV = None
      if unitType == "BQML" and doSUV and factorValid:
        sequenceSUV = suvbwFactor

      if unitType == "BQML" and doSUV and sequenceSUV is None:
        logging.warning("[dPET] SUV conversion disabled: series lacks a valid START/ADMIN SUVbw factor.")

      if unitType == "SUV" and not factorValid:
        logging.warning("[dPET] SUVbw values loaded, but no validated inverse SUVbw factor is available; Bq/mL conversion will be disabled in DynamicPET.")

      for fi in range(nFrames):
        progressbar.value = fi
        slicer.app.processEvents()
        if progressbar.wasCanceled:
          break
        frameFiles = files[fi*filesPerFrame:(fi+1)*filesPerFrame]
        svs = scalarVolumePlugin.examine([frameFiles])
        if not svs:
          raise RuntimeError("Frame parse failed")
        frameNode = scalarVolumePlugin.load(svs[0])
        proxy = slicer.mrmlScene.AddNewNodeByClass(frameNode.GetClassName())
        idx = str(fi)
        volumeSequenceNode.SetDataNodeAtValue(proxy, idx)
        slicer.mrmlScene.RemoveNode(proxy)
        volumeSequenceNode.UpdateDataNodeAtValue(frameNode, idx, True)

        # get stored data node for this index (some Slicer versions return the actual node)
        seqDataNode = volumeSequenceNode.GetDataNodeAtValue(idx) or frameNode
        seqDataNode.SetAttribute("dPETImporter.LoadedBy", "dPETImporterPlugin")

        # DICOMScalarVolumePlugin has now established the final voxel K order
        # and its DICOM.instanceUIDs attribute follows that order exactly.
        # Canonicalize and independently verify the persisted source references
        # against the loaded slice positions before accepting per-slice mapping.
        if (dynamicRTProvenance is not None
            and fi < len(dynamicRTProvenance.get('frames', []))):
          self._canonicalizeLoadedFrameProvenance(
            frameNode, dynamicRTProvenance['frames'][fi])

        # set per-frame attributes on the stored volume node
        if fi < len(frame_durations) and frame_durations[fi] not in (None, ""):
          seqDataNode.SetAttribute('dPET.Duration', str(frame_durations[fi]))
          volumeSequenceNode.SetAttribute(f'dPET.Frame.{idx}.Duration', str(frame_durations[fi]))

        if fi < len(frame_times) and frame_times[fi]:
          seqDataNode.SetAttribute('dPET.AcquisitionTime', str(frame_times[fi]))
          volumeSequenceNode.SetAttribute(f'dPET.Frame.{idx}.AcquisitionTime', str(frame_times[fi]))

        # --- ValueType + SUV handling ---
        if unitType == "BQML":
          if doSUV and (sequenceSUV is not None):
            ok = self._multiplyVolumeByConstant(seqDataNode, sequenceSUV)
            if ok:
              seqDataNode.SetAttribute('dPET.ValueType', 'SUVbw')
              sequenceValueType = "SUVbw"
            else:
              seqDataNode.SetAttribute('dPET.ValueType', 'BQML')
          else:
            seqDataNode.SetAttribute('dPET.ValueType', 'BQML')
        else:
          # DICOM is already SUVbw. Do not rescale the voxel values.
          seqDataNode.SetAttribute('dPET.ValueType', 'SUVbw')

        # Preserve the physical conversion factor even when the stored values
        # remain in BQML or arrived already as SUVbw. Never use 1.0 as a fake
        # inverse conversion factor.
        if factorValid:
          seqDataNode.SetAttribute('dPET.SUVbwFactor', str(suvbwFactor))
          seqDataNode.SetAttribute('dPET.SUVbwFactorValid', '1')
          volumeSequenceNode.SetAttribute(f'dPET.Frame.{idx}.SUVbwFactor', str(suvbwFactor))
          volumeSequenceNode.SetAttribute(f'dPET.Frame.{idx}.SUVbwFactorValid', '1')
        else:
          seqDataNode.SetAttribute('dPET.SUVbwFactorValid', '0')
          volumeSequenceNode.SetAttribute(f'dPET.Frame.{idx}.SUVbwFactorValid', '0')

        # cleanup temporary nodes created by scalarVolumePlugin.load
        if frameNode.GetDisplayNode():
          slicer.mrmlScene.RemoveNode(frameNode.GetDisplayNode())
        if frameNode.GetStorageNode():
          slicer.mrmlScene.RemoveNode(frameNode.GetStorageNode())
        if frameNode is not seqDataNode:
          slicer.mrmlScene.RemoveNode(frameNode)

      # Persist finalized source references only after every loaded frame has
      # been reconciled with Slicer's canonical DICOM.instanceUIDs ordering.
      if dynamicRTProvenance is not None:
        provenanceJson = json.dumps(
          dynamicRTProvenance, separators=(',', ':'), ensure_ascii=False)
        volumeSequenceNode.SetAttribute('dPET.DICOM.FrameReferences', provenanceJson)
        volumeSequenceNode.SetAttribute(
          'dPET.DICOM.ProvenanceSchemaVersion',
          str(dynamicRTProvenance.get('schemaVersion', 2)))

      # create browser and show sequence
      volumeSequenceNode.SetAttribute('dPET.ValueType', sequenceValueType)

      if factorValid:
        volumeSequenceNode.SetAttribute("dPET.SUVbwFactor", str(suvbwFactor))
        volumeSequenceNode.SetAttribute("dPET.SUVbwFactorValid", "1")
      else:
        volumeSequenceNode.SetAttribute("dPET.SUVbwFactorValid", "0")

      browser = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSequenceBrowserNode', slicer.mrmlScene.GenerateUniqueName(baseName + " browser"))
      if slicer.modules.sequences.widgetRepresentation():
        slicer.modules.sequences.widgetRepresentation().setActiveBrowserNode(browser)
      browser.SetAndObserveMasterSequenceNodeID(volumeSequenceNode.GetID())
      browser.SetSaveChanges(volumeSequenceNode, True)
      browser.SetOverwriteProxyName(volumeSequenceNode, True)
      browser.SetAttribute("dPETImporter.LoadedBy", "dPETImporterPlugin")
      browser.SetAttribute("dPETImporter.SequenceNodeID", volumeSequenceNode.GetID())
      browser.SetAttribute('dPET.ValueType', sequenceValueType)
      proxyVol = browser.GetProxyNode(volumeSequenceNode)
      if proxyVol:
        proxyVol.SetAttribute("dPETImporter.LoadedBy", "dPETImporterPlugin")
        proxyVol.SetAttribute('dPET.ValueType', sequenceValueType)
        self._applyKineticMetadata(proxyVol, kineticMetadata)

        # Preserve injection/acquisition provenance on the proxy too.  The
        # sequence already receives all mvNode attributes above, but DynamicPET
        # and other modules frequently inspect the proxy scalar volume directly.
        for attrName in (
          'RadiopharmaceuticalStartDateTime',
          'RadiopharmaceuticalStartTime',
          'RadionuclideStartDateTime',
          'RadionuclideTotalDose',
          'RadionuclideHalfLife',
          'PatientWeight',
          'DecayCorrection',
          'CorrectedImage',
          'dPET.FirstFrameAcquisitionDateTime',
          'dPET.InjectionDateTimeSource',
          'dPET.InjectionToAcquisitionOffsetSec'):
          attrValue = volumeSequenceNode.GetAttribute(attrName)
          if attrValue not in (None, ''):
            proxyVol.SetAttribute(attrName, attrValue)
        if factorValid:
          proxyVol.SetAttribute('dPET.SUVbwFactor', str(suvbwFactor))
          proxyVol.SetAttribute('dPET.SUVbwFactorValid', '1')
        else:
          proxyVol.SetAttribute('dPET.SUVbwFactorValid', '0')
        self._setProxyQuantityUnits(proxyVol, sequenceValueType)

        disp = proxyVol.GetDisplayNode() or proxyVol.CreateDefaultDisplayNodes()
        disp = proxyVol.GetDisplayNode()
        if disp:
          self.setPetDicomLUT(proxyVol)
        appLogic = slicer.app.applicationLogic()
        selNode = appLogic.GetSelectionNode()
        selNode.SetReferenceActiveVolumeID(proxyVol.GetID())
        appLogic.PropagateVolumeSelection()
      # add to subject hierarchy
      self.addSeriesInSubjectHierarchy(loadable, proxyVol if proxyVol else volumeSequenceNode)

      _dpetWatchPETSequenceDisplay(browser, volumeSequenceNode, proxyVol)

      # The proxy relationship is now fully established. Disable automatic
      # proxy renaming only at this final stage (equivalent to unchecking
      # "Rename" manually in the Sequences module after loading). This does
      # not affect playback or synchronization; it only keeps the proxy name
      # stable while the selected item changes.
      if proxyVol:
        browser.SetOverwriteProxyName(volumeSequenceNode, False)
        proxyVol.SetName(baseName)
        try:
          shNode = slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(
              slicer.mrmlScene)
          proxyItemID = shNode.GetItemByDataNode(proxyVol) if shNode else 0
          if (shNode and proxyItemID !=
              slicer.vtkMRMLSubjectHierarchyNode.INVALID_ITEM_ID):
            shNode.SetItemName(proxyItemID, baseName)
        except Exception:
          pass

      return volumeSequenceNode
    except Exception as e:
      logging.error(f"dPET import failed: {e}")
      import traceback
      traceback.print_exc()
      return None
    finally:
      progressbar.close()


# -----------------------------------------------------------------------------
# SlicerDynamicPET DICOM Parametric Map importer
# -----------------------------------------------------------------------------

class dPETParametricMapPluginClass(DICOMPlugin):

  PARAMETRIC_MAP_SOP_CLASS_UID = \
    "1.2.840.10008.5.1.4.1.1.30"

  def __init__(self):
    super().__init__()

    self.loadType = "SlicerDynamicPET Parametric Map"

    self.tags['sopClassUID'] = "0008,0016"
    self.tags['sopInstanceUID'] = "0008,0018"
    self.tags['seriesInstanceUID'] = "0020,000E"
    self.tags['studyInstanceUID'] = "0020,000D"
    self.tags['seriesDescription'] = "0008,103E"
    self.tags['seriesNumber'] = "0020,0011"
    self.tags['modality'] = "0008,0060"


  def _ensureHighdicom(self):
    try:
      import highdicom as hd
      return hd

    except ImportError:
      logging.info(
        "[dPET PM] Installing highdicom without dependencies..."
      )

      # Slicer has a shared Python environment. Do not let pip upgrade or
      # replace Slicer's pydicom/numpy/etc. while enabling PM import.
      slicer.util.pip_install(
        "--no-deps highdicom"
      )

      import importlib
      importlib.invalidate_caches()

      import highdicom as hd
      return hd


  def examine(self, fileLists):
    """
    Offer DICOM Parametric Map Storage objects using our
    dedicated PM loader instead of Slicer's generic scalar
    volume loader.
    """

    loadables = []

    for files in fileLists:
      if not files:
        continue

      for filePath in files:

        sopClassUID = (
          slicer.dicomDatabase.fileValue(
            filePath,
            self.tags['sopClassUID']
          )
          or ""
        )

        if (
          sopClassUID
          != self.PARAMETRIC_MAP_SOP_CLASS_UID
        ):
          continue

        seriesDescription = (
          slicer.dicomDatabase.fileValue(
            filePath,
            self.tags['seriesDescription']
          )
          or
          "SlicerDynamicPET Parametric Map"
        )

        seriesNumber = (
          slicer.dicomDatabase.fileValue(
            filePath,
            self.tags['seriesNumber']
          )
          or ""
        )

        loadable = DICOMLoadable()

        # One PM SOP instance is one loadable.
        loadable.files = [filePath]

        if seriesNumber:
          loadable.name = (
            f"{seriesNumber}: "
            f"{seriesDescription}"
          )
        else:
          loadable.name = seriesDescription

        loadable.tooltip = (
          "DICOM Parametric Map loaded by "
          "SlicerDynamicPET/dPETImporter"
        )

        # Give this dedicated PM loader strong preference
        # over the generic Scalar Volume DICOM plugin.
        loadable.selected = True
        loadable.confidence = 1.0

        loadables.append(loadable)

    return loadables


  def _firstRealWorldValueMapping(self, pm):
    """
    Return the first Real World Value Mapping item.

    PM commonly stores this in the Shared Functional Groups,
    but handle per-frame storage too.
    """

    shared = getattr(
      pm,
      "SharedFunctionalGroupsSequence",
      None
    )

    if shared and len(shared) > 0:
      mappingSequence = getattr(
        shared[0],
        "RealWorldValueMappingSequence",
        None
      )

      if (
        mappingSequence is not None
        and len(mappingSequence) > 0
      ):
        return mappingSequence[0]

    perFrame = getattr(
      pm,
      "PerFrameFunctionalGroupsSequence",
      None
    )

    if perFrame:
      for functionalGroup in perFrame:

        mappingSequence = getattr(
          functionalGroup,
          "RealWorldValueMappingSequence",
          None
        )

        if (
          mappingSequence is not None
          and len(mappingSequence) > 0
        ):
          return mappingSequence[0]

    return None


  def _makeCodedEntry(
      self,
      codeValue,
      codingScheme,
      codeMeaning):

    CodedEntryClass = (
      getattr(vtk, "vtkCodedEntry", None)
      or
      getattr(slicer, "vtkCodedEntry", None)
    )

    if CodedEntryClass is None:
      return None

    entry = CodedEntryClass()

    entry.SetValueSchemeMeaning(
      str(codeValue),
      str(codingScheme),
      str(codeMeaning)
    )

    return entry


  def _setParametricMapQuantityMetadata(
      self,
      volumeNode,
      pm):

    mapping = self._firstRealWorldValueMapping(
      pm
    )

    if mapping is None:
      return

    quantityCode = str(
      getattr(
        mapping,
        "LUTLabel",
        ""
      )
      or ""
    )

    quantityMeaning = str(
      getattr(
        mapping,
        "LUTExplanation",
        quantityCode
      )
      or quantityCode
    )

    unitCode = ""
    unitScheme = ""
    unitMeaning = ""

    unitsSequence = getattr(
      mapping,
      "MeasurementUnitsCodeSequence",
      None
    )

    if (
      unitsSequence is not None
      and len(unitsSequence) > 0
    ):
      units = unitsSequence[0]

      unitCode = str(
        getattr(
          units,
          "CodeValue",
          ""
        )
        or ""
      )

      unitScheme = str(
        getattr(
          units,
          "CodingSchemeDesignator",
          ""
        )
        or ""
      )

      unitMeaning = str(
        getattr(
          units,
          "CodeMeaning",
          unitCode
        )
        or unitCode
      )

    # Keep textual metadata even on Slicer builds where
    # vtkCodedEntry is unavailable.
    if quantityCode:
      volumeNode.SetAttribute(
        "SlicerDynamicPET.QuantityCode",
        quantityCode
      )

    if quantityMeaning:
      volumeNode.SetAttribute(
        "SlicerDynamicPET.QuantityMeaning",
        quantityMeaning
      )

    if unitCode:
      volumeNode.SetAttribute(
        "SlicerDynamicPET.UnitCode",
        unitCode
      )

    if unitScheme:
      volumeNode.SetAttribute(
        "SlicerDynamicPET.UnitCodingScheme",
        unitScheme
      )

    if unitMeaning:
      volumeNode.SetAttribute(
        "SlicerDynamicPET.UnitMeaning",
        unitMeaning
      )

    # Also populate Slicer's quantitative voxel metadata.
    if (
      quantityCode
      and
      hasattr(
        volumeNode,
        "SetVoxelValueQuantity"
      )
    ):
      quantity = self._makeCodedEntry(
        quantityCode,
        "99SDPET",
        quantityMeaning
      )

      if quantity is not None:
        volumeNode.SetVoxelValueQuantity(
          quantity
        )

    if (
      unitCode
      and
      hasattr(
        volumeNode,
        "SetVoxelValueUnits"
      )
    ):
      units = self._makeCodedEntry(
        unitCode,
        unitScheme or "UCUM",
        unitMeaning
      )

      if units is not None:
        volumeNode.SetVoxelValueUnits(
          units
        )


  @staticmethod
  def _firstNumericValue(value):
    if value is None:
      return None
    try:
      if isinstance(value, (list, tuple)):
        if len(value) == 0:
          return None
        return float(value[0])
      # pydicom MultiValue behaves as an iterable but should not be confused
      # with a numeric string.
      if not isinstance(value, (str, bytes)) and hasattr(value, "__len__"):
        if len(value) == 0:
          return None
        return float(value[0])
      text = str(value).strip()
      if "\\" in text:
        text = text.split("\\", 1)[0]
      return float(text)
    except Exception:
      return None


  def _parametricMapStoredWindow(self, pm):
    """Return (minimum, maximum) from DICOM VOI window, if present."""
    import numpy as np

    def fromDataset(dataset):
      if dataset is None:
        return None
      center = self._firstNumericValue(
        getattr(dataset, "WindowCenter", None)
      )
      width = self._firstNumericValue(
        getattr(dataset, "WindowWidth", None)
      )
      if center is None or width is None or width <= 0.0:
        return None
      minimum = center - 0.5 * width
      maximum = center + 0.5 * width
      if not (np.isfinite(minimum) and np.isfinite(maximum)):
        return None
      if maximum <= minimum:
        return None
      return float(minimum), float(maximum)

    result = fromDataset(pm)
    if result is not None:
      return result

    shared = getattr(pm, "SharedFunctionalGroupsSequence", None)
    if shared:
      voi = getattr(shared[0], "FrameVOILUTSequence", None)
      if voi:
        result = fromDataset(voi[0])
        if result is not None:
          return result

    perFrame = getattr(pm, "PerFrameFunctionalGroupsSequence", None)
    if perFrame:
      voi = getattr(perFrame[0], "FrameVOILUTSequence", None)
      if voi:
        result = fromDataset(voi[0])
        if result is not None:
          return result

    return None


  @staticmethod
  def _fallbackParametricDisplayRange(pixelArray, quantityCode):
    """Robust fallback only for PMs that do not carry usable DICOM W/L."""
    import numpy as np
    values = np.asarray(pixelArray, dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size == 0:
      return None

    if quantityCode == "SDP_VB":
      return 0.0, 1.0, "Physical[0,1]"

    q25, q75 = np.percentile(values, [25.0, 75.0])
    centralWidth = max(0.0, float(q75 - q25))

    if (
      quantityCode == "SDP_TD"
      or (
        quantityCode == "SDP_INT"
        and float(np.min(values)) < 0.0
        and float(np.max(values)) > 0.0
      )
    ):
      absValues = np.abs(values)
      selected = None
      for percentile, label in (
          (99.0, "SymmetricP99Abs"),
          (98.0, "SymmetricP98Abs"),
          (95.0, "SymmetricP95Abs")):
        halfRange = float(np.percentile(absValues, percentile))
        if not np.isfinite(halfRange) or halfRange <= 0.0:
          continue
        selected = (-halfRange, halfRange, label)
        occupancy = centralWidth / max(2.0 * halfRange, 1.0e-12)
        if occupancy >= 0.20:
          break
      return selected

    for lowPercentile, highPercentile, label in (
        (1.0, 99.0, "P1-P99"),
        (2.0, 98.0, "P2-P98"),
        (5.0, 95.0, "P5-P95")):
      minimum, maximum = np.percentile(
        values,
        [lowPercentile, highPercentile]
      )
      minimum = float(minimum)
      maximum = float(maximum)
      if not (
        np.isfinite(minimum)
        and np.isfinite(maximum)
        and maximum > minimum
      ):
        continue
      occupancy = centralWidth / max(maximum - minimum, 1.0e-12)
      selected = (minimum, maximum, label)
      if occupancy >= 0.20 or label == "P5-P95":
        return selected

    minimum = float(np.min(values))
    maximum = float(np.max(values))
    if maximum <= minimum:
      scale = max(1.0, abs(minimum), abs(maximum))
      epsilon = 1.0e-6 * scale
      minimum -= 0.5 * epsilon
      maximum += 0.5 * epsilon
    return minimum, maximum, "FiniteMinMaxFallback"


  def _applyParametricMapDisplayPreset(
      self,
      volumeNode,
      pm,
      pixelArray):
    """Restore SlicerDynamicPET display intent without changing voxel values."""
    import numpy as np
    if volumeNode is None:
      return

    volumeNode.CreateDefaultDisplayNodes()
    displayNode = volumeNode.GetDisplayNode()
    if displayNode is None:
      return

    quantityCode = str(
      volumeNode.GetAttribute("SlicerDynamicPET.QuantityCode")
      or ""
    )

    finiteValues = np.asarray(pixelArray, dtype=np.float64)
    finiteValues = finiteValues[np.isfinite(finiteValues)]
    interceptSpansZero = (
      quantityCode == "SDP_INT"
      and finiteValues.size > 0
      and float(np.min(finiteValues)) < 0.0
      and float(np.max(finiteValues)) > 0.0
    )

    if quantityCode in (
        "SDP_AIC",
        "SDP_BIC",
        "SDP_MASE",
        "SDP_CHI2"):
      colorNodeID = "vtkMRMLColorTableNodeFileMagma.txt"
    elif quantityCode == "SDP_TD" or interceptSpansZero:
      colorNodeID = "vtkMRMLColorTableNodeFileDivergingBlueRed.txt"
    else:
      colorNodeID = "vtkMRMLColorTableNodeFileViridis.txt"

    # vb is intentionally always displayed over its physical [0,1] domain.
    if quantityCode == "SDP_VB":
      displayRange = (0.0, 1.0)
      policy = "Physical[0,1]"
    else:
      displayRange = self._parametricMapStoredWindow(pm)
      policy = "DICOMWindowLevel"

      if displayRange is None:
        fallback = self._fallbackParametricDisplayRange(
          pixelArray,
          quantityCode
        )
        if fallback is None:
          return
        displayRange = fallback[:2]
        policy = fallback[2]

    minimum, maximum = displayRange
    if not (
      np.isfinite(minimum)
      and np.isfinite(maximum)
      and maximum > minimum
    ):
      return

    displayNode.SetAutoWindowLevel(False)
    displayNode.SetWindowLevelMinMax(
      float(minimum),
      float(maximum)
    )
    displayNode.SetInterpolate(True)
    displayNode.SetInvertDisplayScalarRange(False)
    displayNode.SetWindowMappingMethod(
      slicer.vtkMRMLScalarVolumeDisplayNode.GetWindowMappingMethodFromString(
        "Linear"
      )
    )
    displayNode.SetApplyThreshold(False)
    displayNode.SetAutoThreshold(False)
    displayNode.SetAndObserveColorNodeID(colorNodeID)

    volumeNode.SetAttribute(
      "SlicerDynamicPET.Display.RangePolicy",
      str(policy)
    )
    volumeNode.SetAttribute(
      "SlicerDynamicPET.Display.WindowMin",
      format(float(minimum), ".16g")
    )
    volumeNode.SetAttribute(
      "SlicerDynamicPET.Display.WindowMax",
      format(float(maximum), ".16g")
    )
    volumeNode.SetAttribute(
      "SlicerDynamicPET.Display.ColorNodeID",
      colorNodeID
    )
    volumeNode.SetAttribute(
      "SlicerDynamicPET.Display.Mapping",
      "Linear"
    )


  def load(self, loadable):
    """
    Load a DICOM Parametric Map as a vtkMRMLScalarVolumeNode.

    Pixel values are returned in real-world units and the
    DICOM LPS geometry is converted explicitly to Slicer RAS.
    """

    volumeNode = None

    try:
      import numpy as np

      hd = self._ensureHighdicom()

      if not loadable.files:
        raise RuntimeError(
          "Parametric Map loadable has no DICOM file."
        )

      filePath = loadable.files[0]

      # --------------------------------------------------------------
      # Read PM using highdicom's dedicated Parametric Map reader.
      # --------------------------------------------------------------

      pmread = getattr(
        getattr(hd, "pm", None),
        "pmread",
        None
      )

      pm = (
        pmread(filePath)
        if callable(pmread)
        else hd.imread(filePath)
      )

      if (
        str(pm.SOPClassUID)
        != self.PARAMETRIC_MAP_SOP_CLASS_UID
      ):
        raise RuntimeError(
          "Selected object is not a DICOM "
          "Parametric Map Storage instance."
        )

      # --------------------------------------------------------------
      # Recover real-world voxel values and spatial geometry.
      #
      # Float Pixel Data intentionally has no PixelRepresentation in
      # DICOM. Some highdicom releases nevertheless access that attribute
      # in Image.get_volume(), so use pydicom's float-aware pixel decoder
      # and highdicom's geometry-only API for floating-point PMs.
      # Integer PMs continue through highdicom.get_volume().
      # --------------------------------------------------------------

      hasFloatPixelData = (
        hasattr(pm, "FloatPixelData")
        or hasattr(pm, "DoubleFloatPixelData")
      )

      if hasFloatPixelData and not hasattr(pm, "PixelRepresentation"):
        storedArray = np.asarray(pm.pixel_array)

        if storedArray.ndim == 2:
          storedArray = storedArray[np.newaxis, :, :]

        if storedArray.ndim != 3:
          raise RuntimeError(
            "Expected a 3D scalar Parametric Map, "
            f"but got stored pixel shape {storedArray.shape}."
          )

        # Apply the PM's real-world value mapping. SlicerDynamicPET writes
        # a linear mapping, but support a standard LUT mapping as well.
        mapping = self._firstRealWorldValueMapping(pm)
        if mapping is None:
          raise RuntimeError(
            "Parametric Map does not contain a "
            "Real World Value Mapping."
          )

        slope = getattr(mapping, "RealWorldValueSlope", None)
        intercept = getattr(mapping, "RealWorldValueIntercept", None)

        if slope is not None and intercept is not None:
          pixelArrayStoredOrder = (
            storedArray.astype(np.float32, copy=False)
            * float(slope)
            + float(intercept)
          )
        else:
          lutData = getattr(mapping, "RealWorldValueLUTData", None)
          firstMapped = getattr(
            mapping,
            "RealWorldValueFirstValueMapped",
            None
          )
          if firstMapped is None:
            firstMapped = getattr(
              mapping,
              "DoubleFloatRealWorldValueFirstValueMapped",
              None
            )

          if lutData is None or firstMapped is None:
            raise RuntimeError(
              "Unsupported Real World Value Mapping in "
              "floating-point Parametric Map."
            )

          if not np.issubdtype(storedArray.dtype, np.integer):
            raise RuntimeError(
              "A LUT-based Real World Value Mapping requires "
              "integer stored pixel values."
            )

          lut = np.asarray(lutData, dtype=np.float32)
          lutIndices = storedArray.astype(np.int64) - int(firstMapped)
          if (
            np.any(lutIndices < 0)
            or np.any(lutIndices >= lut.size)
          ):
            raise RuntimeError(
              "Stored Parametric Map pixels fall outside the "
              "Real World Value LUT range."
            )
          pixelArrayStoredOrder = lut[lutIndices]

        # Geometry can be obtained without decoding/transformation of pixels.
        geometry = pm.get_volume_geometry()
        if geometry is None:
          raise RuntimeError(
            "Parametric Map frames do not form a regular 3D volume."
          )

        kjiToLPS = np.asarray(
          geometry.affine,
          dtype=np.float64
        )

        if kjiToLPS.shape != (4, 4):
          raise RuntimeError(
            "Invalid Parametric Map affine matrix."
          )

        # pm.pixel_array is in stored DICOM frame order, while the volume
        # geometry is in spatial K order. Match each spatial plane back to
        # its stored frame by Image Position (Patient) rather than assuming
        # those orders are identical.
        perFrame = getattr(
          pm,
          "PerFrameFunctionalGroupsSequence",
          None
        )
        if perFrame is None or len(perFrame) != storedArray.shape[0]:
          raise RuntimeError(
            "Parametric Map frame geometry is incomplete."
          )

        storedPositions = []
        for frameGroup in perFrame:
          positionSequence = getattr(
            frameGroup,
            "PlanePositionSequence",
            None
          )
          if not positionSequence:
            raise RuntimeError(
              "Parametric Map frame is missing Plane Position."
            )
          storedPositions.append(
            np.asarray(
              positionSequence[0].ImagePositionPatient,
              dtype=np.float64
            )
          )

        spatialPositions = []
        for positionSequence in geometry.get_plane_positions():
          spatialPositions.append(
            np.asarray(
              positionSequence[0].ImagePositionPatient,
              dtype=np.float64
            )
          )

        if len(spatialPositions) != len(storedPositions):
          raise RuntimeError(
            "Parametric Map spatial frame count does not match "
            "the stored pixel frame count."
          )

        frameOrder = []
        usedStoredFrames = set()
        for spatialPosition in spatialPositions:
          distances = [
            float(np.linalg.norm(position - spatialPosition))
            for position in storedPositions
          ]
          storedIndex = int(np.argmin(distances))
          if (
            distances[storedIndex] > 1.0e-3
            or storedIndex in usedStoredFrames
          ):
            raise RuntimeError(
              "Could not match Parametric Map stored frames to "
              "their spatial volume positions."
            )
          usedStoredFrames.add(storedIndex)
          frameOrder.append(storedIndex)

        pixelArray = np.asarray(
          pixelArrayStoredOrder[frameOrder, :, :],
          dtype=np.float32
        )

      else:
        parametricVolume = pm.get_volume(
          dtype=np.float32,
          apply_real_world_transform=True
        )

        pixelArray = np.asarray(
          parametricVolume.array,
          dtype=np.float32
        )

        if pixelArray.ndim != 3:
          raise RuntimeError(
            "Expected a 3D scalar Parametric Map, "
            f"but got shape {pixelArray.shape}."
          )

        kjiToLPS = np.asarray(
          parametricVolume.affine,
          dtype=np.float64
        )

        if kjiToLPS.shape != (4, 4):
          raise RuntimeError(
            "Invalid Parametric Map affine matrix."
          )

      # --------------------------------------------------------------
      # Convert highdicom KJI indexing to Slicer IJK indexing.
      #
      # highdicom:
      #   axis 0 -> K
      #   axis 1 -> J
      #   axis 2 -> I
      #
      # Slicer:
      #   matrix column 0 -> I
      #   matrix column 1 -> J
      #   matrix column 2 -> K
      # --------------------------------------------------------------

      ijkToLPS = np.eye(
        4,
        dtype=np.float64
      )

      ijkToLPS[:3, 0] = \
        kjiToLPS[:3, 2]

      ijkToLPS[:3, 1] = \
        kjiToLPS[:3, 1]

      ijkToLPS[:3, 2] = \
        kjiToLPS[:3, 0]

      ijkToLPS[:3, 3] = \
        kjiToLPS[:3, 3]

      # DICOM patient coordinates are LPS.
      # Slicer uses RAS.
      lpsToRAS = np.array(
        [
          [-1.0,  0.0, 0.0, 0.0],
          [ 0.0, -1.0, 0.0, 0.0],
          [ 0.0,  0.0, 1.0, 0.0],
          [ 0.0,  0.0, 0.0, 1.0],
        ],
        dtype=np.float64
      )

      ijkToRASArray = (
        lpsToRAS @ ijkToLPS
      )

      # --------------------------------------------------------------
      # Create Slicer scalar volume
      # --------------------------------------------------------------

      seriesDescription = str(
        getattr(
          pm,
          "SeriesDescription",
          ""
        )
        or
        "SlicerDynamicPET Parametric Map"
      )

      nodeName = (
        slicer.mrmlScene.GenerateUniqueName(
          seriesDescription
        )
      )

      volumeNode = (
        slicer.mrmlScene.AddNewNodeByClass(
          "vtkMRMLScalarVolumeNode",
          nodeName
        )
      )

      slicer.util.updateVolumeFromArray(
        volumeNode,
        pixelArray
      )

      ijkToRAS = vtk.vtkMatrix4x4()

      for row in range(4):
        for column in range(4):
          ijkToRAS.SetElement(
            row,
            column,
            float(
              ijkToRASArray[
                row,
                column
              ]
            )
          )

      volumeNode.SetIJKToRASMatrix(
        ijkToRAS
      )

      # Match the protection used by Slicer's standard
      # scalar-volume DICOM loader: keep IJK right-handed
      # while preserving physical geometry.
      if not slicer.vtkMRMLVolumeNode.IsIJKCoordinateSystemRightHanded(
          ijkToRAS):

        slicer.vtkMRMLVolumeNode.ReverseSliceOrder(
          volumeNode.GetImageData(),
          ijkToRAS
        )

        volumeNode.SetIJKToRASMatrix(
          ijkToRAS
        )

      # --------------------------------------------------------------
      # DICOM provenance
      # --------------------------------------------------------------

      volumeNode.SetAttribute(
        "DICOM.instanceUIDs",
        str(pm.SOPInstanceUID)
      )

      volumeNode.SetAttribute(
        "DICOM.SeriesInstanceUID",
        str(pm.SeriesInstanceUID)
      )

      volumeNode.SetAttribute(
        "DICOM.StudyInstanceUID",
        str(pm.StudyInstanceUID)
      )

      volumeNode.SetAttribute(
        "DICOM.SOPClassUID",
        str(pm.SOPClassUID)
      )

      if hasattr(
          pm,
          "FrameOfReferenceUID"):
        volumeNode.SetAttribute(
          "DICOM.FrameOfReferenceUID",
          str(pm.FrameOfReferenceUID)
        )

      volumeNode.SetAttribute(
        "SlicerDynamicPET.ResultType",
        "ParametricMap"
      )

      volumeNode.SetAttribute(
        "SlicerDynamicPET.ImportedFromDICOM",
        "1"
      )

      derivationDescription = str(
        getattr(
          pm,
          "DerivationDescription",
          ""
        )
        or ""
      )

      if derivationDescription:
        volumeNode.SetAttribute(
          "SlicerDynamicPET.DerivationDescription",
          derivationDescription
        )

      # Quantity and units from the PM Real World
      # Value Mapping sequence.
      self._setParametricMapQuantityMetadata(
        volumeNode,
        pm
      )

      # --------------------------------------------------------------
      # Display
      #
      # Prefer the robust window/level persisted in the DICOM PM.  The
      # quantity code selects a Slicer LUT, while mapping remains linear and
      # thresholding remains disabled.
      # --------------------------------------------------------------

      self._applyParametricMapDisplayPreset(
        volumeNode,
        pm,
        pixelArray
      )

      # --------------------------------------------------------------
      # Put it into the DICOM patient/study hierarchy.
      #
      # This uses the PM object's own patient/study DICOM
      # metadata and therefore keeps it under the same study.
      # --------------------------------------------------------------

      self.addSeriesInSubjectHierarchy(
        loadable,
        volumeNode
      )

      # Make loaded PM visible as active volume.
      appLogic = (
        slicer.app.applicationLogic()
      )

      if appLogic:
        selectionNode = (
          appLogic.GetSelectionNode()
        )

        selectionNode.SetReferenceActiveVolumeID(
          volumeNode.GetID()
        )

        appLogic.PropagateVolumeSelection()

      logging.info(
        "[dPET PM] Loaded Parametric Map "
        f"'{seriesDescription}', "
        f"shape={pixelArray.shape}, "
        f"SOPInstanceUID={pm.SOPInstanceUID}"
      )

      return volumeNode

    except Exception as e:

      logging.error(
        f"[dPET PM] Parametric Map load failed: {e}"
      )

      import traceback
      traceback.print_exc()

      if volumeNode is not None:
        try:
          slicer.mrmlScene.RemoveNode(
            volumeNode
          )
        except Exception:
          pass

      return None


#
# dPETImporterPlugin
#

class dPETImporterPlugin:
  """
  This class is the 'hook' for slicer to detect and recognize the plugin
  as a loadable scripted module
  """
  def __init__(self, parent):
    parent.title = "DICOM dPET Import Plugin"
    parent.categories = ["Developer Tools.DICOM Plugins"]
    parent.contributors = ["Daniele Dall'Olio, University of Bologna"]
    parent.helpText = """
    Plugin to the DICOM Module to parse and load dPET data from DICOM files.
    No module interface here, only in the DICOM module
    """
    parent.acknowledgementText = """
    This DICOM Plugin was developed by
    Daniele Dall'Olio, University of Bologna.
    """

    # don't show this module - it only appears in the DICOM module
    parent.hidden = True

    # Add this extension to the DICOM module's list for discovery when the module
    # is created.  Since this module may be discovered before DICOM itself,
    # create the list if it doesn't already exist.
    try:
      slicer.modules.dicomPlugins
    except AttributeError:
      slicer.modules.dicomPlugins = {}

    slicer.modules.dicomPlugins[
      'dPETImporterPlugin'
    ] = dPETImporterPluginClass

    slicer.modules.dicomPlugins[
      'dPETParametricMapPlugin'
    ] = dPETParametricMapPluginClass

#
#

class dPETImporterPluginWidget:
  def __init__(self, parent = None):
    self.parent = parent

  def setup(self):
    # don't display anything for this widget - it will be hidden anyway
    pass

  def enter(self):
    pass

  def exit(self):
    pass
