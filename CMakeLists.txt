cmake_minimum_required(VERSION 3.16.3)
set(EXTENSION_NAME dPETImporter)

if(NOT Slicer_SOURCE_DIR)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
  set(EXTENSION_HOMEPAGE "https://github.com/DanieleDallOlio/dPETImporter")
  set(EXTENSION_CATEGORY "DICOM")
  set(EXTENSION_STATUS "Beta")
  set(EXTENSION_DESCRIPTION "DICOM plugin to import dynamic PET as Volume Sequence")
  set(EXTENSION_DEPENDS NA)
  include(SlicerEnableExtensionTesting)
endif()

add_subdirectory(dPETImporter)

if(NOT Slicer_SOURCE_DIR)
  include(${Slicer_EXTENSION_CPACK})
endif()
