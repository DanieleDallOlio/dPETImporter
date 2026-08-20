cmake_minimum_required(VERSION 3.16.3)

project(dPETImporter)

#-----------------------------------------------------------------------------
# Extension meta-information
set(
  EXTENSION_HOMEPAGE
  "https://github.com/DanieleDallOlio/dPETImporter"
)

set(
  EXTENSION_CATEGORY
  "DICOM"
)

set(
  EXTENSION_CONTRIBUTORS
  "Daniele Dall'Olio (University of Bologna)"
)

set(
  EXTENSION_DESCRIPTION
  "DICOM plugin for 3D Slicer that recognizes dynamic PET studies and loads them as time-resolved volume sequences."
)

set(
  EXTENSION_ICONURL
  "https://raw.githubusercontent.com/DanieleDallOlio/dPETImporter/main/dPETImporter.png"
)

set(
  EXTENSION_SCREENSHOTURLS
  "https://raw.githubusercontent.com/DanieleDallOlio/dPETImporter/main/Screenshots/dPETImporter-DICOM.png"
)

set(
  EXTENSION_DEPENDS
  "NA"
)

#-----------------------------------------------------------------------------
# Find Slicer
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(dPETImporter)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
