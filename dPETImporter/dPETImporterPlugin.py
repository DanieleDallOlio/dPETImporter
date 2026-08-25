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


def compute_suvbw_for_start(mvNode):
  """
  Compute a single SUVbw factor for the whole dynamic sequence.

  Enforced convention:
    - DecayCorrection must be START
    - START is taken as the acquisition datetime of the first frame

  Returns float or None.
  """
  try:
    rstart = mvNode.GetAttribute('RadionuclideStartDateTime')
    first_frame_dt = mvNode.GetAttribute('dPET.FirstFrameAcquisitionDateTime')

    half_life = float(mvNode.GetAttribute('RadionuclideHalfLife'))
    weight = float(mvNode.GetAttribute('PatientWeight'))
    totalDose = float(mvNode.GetAttribute('RadionuclideTotalDose'))

    decayCorrection = (mvNode.GetAttribute('DecayCorrection') or '').upper()
    correctedImage = (mvNode.GetAttribute('CorrectedImage') or mvNode.GetAttribute('correctedImage') or '').upper()

    if not (rstart and first_frame_dt and half_life > 0 and weight > 0 and totalDose > 0):
      return None

    # Require START
    if decayCorrection != "START":
      return None

    # Conservative check for corrected image
    if not (("ATTN" in correctedImage) and ("DECY" in correctedImage or "DECAY" in correctedImage)):
      return None

    start_dt = _parse_dicom_datetime_value(rstart)
    ref_dt = _parse_dicom_datetime_value(first_frame_dt)
    if start_dt is None or ref_dt is None:
      return None

    decay_seconds = (ref_dt - start_dt).total_seconds()

    dose_kbq = totalDose * 0.001
    if dose_kbq <= 0:
      return None

    decayedDose = dose_kbq * 2 ** (-decay_seconds / half_life)
    if decayedDose <= 0:
      return None

    suvbw = weight / decayedDose
    return float(suvbw)
  except Exception:
    return None



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


  def examine(self,fileLists):
    if not settingsValue("DICOM/dPETImporterEnabled", True, converter=toBool):
      return []

    self.detailedLogging = settingsValue('DICOM/detailedLogging', False, converter=toBool)
    loadables = []
    allfiles = []
    for files in fileLists:
      loadables += self.examineFiles(files)
      allfiles += files

    if (not loadables) and len(allfiles)>len(files):
      loadables += self.examineFilesMultiseries(allfiles)

    # --- annotate each loadable (MV and Sequence) with 2D/3D + SUV parsing status ---
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
      if len(missing) == 0:
        suvTxt = "SUV:ok"
      else:
        # keep it short in tooltip; full list in debug log if you want
        suvTxt = "SUV:missing"

      suffix = f" [{ft}; {suvTxt}]"

      # name: short, tooltip: longer
      if loadable.name:
        loadable.name = loadable.name + suffix
      else:
        loadable.name = suffix.strip()

      if loadable.tooltip:
        if len(missing) == 0:
          loadable.tooltip = loadable.tooltip + suffix
        else:
          loadable.tooltip = loadable.tooltip + suffix + f" (missing: {', '.join(missing)})"
      else:
        loadable.tooltip = suffix.strip()

    if hasattr(slicer.modules, 'sequences'):
      seqLoadables = []
      for loadable in loadables:
        seqL = DICOMLoadable()
        seqL.files = loadable.files
        seqL.multivolume = loadable.multivolume
        seqL.selected = loadable.selected
        seqL.confidence = loadable.confidence
        seqL.loadAsVolumeSequence = True
        # set name & tooltip (prefer the sequence label)
        seqL.tooltip = (loadable.tooltip or '').replace('MultiVolume', 'Volume Sequence')
        seqL.name = (loadable.name or '').replace('MultiVolume', 'Volume Sequence')
        annotate(seqL, isSequence=True)
        seqLoadables.append(seqL)
      # loadables[0:0] = seqLoadables
      return seqLoadables
    return []

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
    if not volumeNode:
      return False

    # Ensure display node exists
    displayNode = volumeNode.GetDisplayNode()
    if displayNode is None:
      volumeNode.CreateDefaultDisplayNodes()
      displayNode = volumeNode.GetDisplayNode()
    if displayNode is None:
      return False

    # Robustly get PET-DICOM color node
    petColorNode = slicer.mrmlScene.GetFirstNodeByName("PET-DICOM")
    if petColorNode is None:
      # If name lookup fails, get the singleton ID for the PETDICOM procedural node
      petColorNodeID = vtk.vtkMRMLColorLogic.GetPETColorNodeID(
        vtk.vtkMRMLPETProceduralColorNode.PETDICOM
      )
      petColorNode = slicer.mrmlScene.GetNodeByID(petColorNodeID)

    if petColorNode is None:
      logging.error("[dPET] PET-DICOM color node not found in scene.")
      return False

    displayNode.SetAndObserveColorNodeID(petColorNode.GetID())
    displayNode.Modified()
    volumeNode.Modified()
    return True

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

  def load(self,loadable):
    """Load as Volume Sequence (only sequence is supported in this minimal plugin)"""
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
      suvbwFactor = compute_suvbw_for_start(mvNode)
      factorValid = suvbwFactor is not None and suvbwFactor > 0

      sequenceSUV = None
      if unitType == "BQML" and doSUV and factorValid:
        sequenceSUV = suvbwFactor

      if unitType == "BQML" and doSUV and sequenceSUV is None:
        logging.warning("[dPET] SUV conversion disabled: series is not usable as decay-corrected-to-first-frame START.")

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

        # Preserve injection/acquisition provenance on the proxy too.  The
        # sequence already receives all mvNode attributes above, but DynamicPET
        # and other modules frequently inspect the proxy scalar volume directly.
        for attrName in (
          'RadiopharmaceuticalStartDateTime',
          'RadiopharmaceuticalStartTime',
          'RadionuclideStartDateTime',
          'RadionuclideTotalDose',
          'RadionuclideHalfLife',
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
          disp.SetAutoWindowLevel(True)
          disp.SetInterpolate(True)
        appLogic = slicer.app.applicationLogic()
        selNode = appLogic.GetSelectionNode()
        selNode.SetReferenceActiveVolumeID(proxyVol.GetID())
        appLogic.PropagateVolumeSelection()
      # add to subject hierarchy
      self.addSeriesInSubjectHierarchy(loadable, proxyVol if proxyVol else volumeSequenceNode)
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
        "[dPET PM] Installing highdicom 0.28.1..."
      )

      slicer.util.pip_install(
        "highdicom==0.28.1"
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
      # Read PM using highdicom
      # --------------------------------------------------------------

      pm = hd.imread(
          filePath
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
      # Recover real-world voxel values and spatial geometry
      #
      # highdicom Volume uses:
      #   array axes = [K, J, I]
      #   affine     = KJI -> DICOM LPS
      # --------------------------------------------------------------

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
      # --------------------------------------------------------------

      volumeNode.CreateDefaultDisplayNodes()

      displayNode = (
        volumeNode.GetDisplayNode()
      )

      if displayNode:
        displayNode.SetAutoWindowLevel(
          True
        )

        displayNode.SetInterpolate(
          True
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
