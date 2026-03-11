import slicer

class dPETImporter:
  def __init__(self, parent):
    parent.title = "DICOM dPET Import Plugin"
    parent.categories = ["Developer Tools.DICOM Plugins"]
    parent.contributors = ["Daniele Dall'Olio (University of Bologna)"]
    parent.helpText = "DICOM plugin to import dynamic PET datasets as volume sequences."
    parent.hidden = True

    from dPETImporterPlugin import dPETImporterPluginClass

    try:
      slicer.modules.dicomPlugins
    except AttributeError:
      slicer.modules.dicomPlugins = {}

    slicer.modules.dicomPlugins['dPETImporterPlugin'] = dPETImporterPluginClass
