import slicer

class dPETImporter:
  def __init__(self, parent):
    parent.title = "DICOM dPET Import Plugin"
    parent.categories = ["Developer Tools.DICOM Plugins"]
    parent.hidden = True

    from dPETImporterPlugin import dPETImporterPluginClass

    try:
      slicer.modules.dicomPlugins
    except AttributeError:
      slicer.modules.dicomPlugins = {}
    slicer.modules.dicomPlugins['dPETImporterPlugin'] = dPETImporterPluginClass
