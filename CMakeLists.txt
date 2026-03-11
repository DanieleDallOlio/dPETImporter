cmake_minimum_required(VERSION 3.16.3)

project(dPETImporter)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://github.com/DanieleDallOlio/dPETImporter")
set(EXTENSION_CATEGORY "DICOM")
set(EXTENSION_CONTRIBUTORS "Daniele Dall'Olio (University of Bologna)")
set(EXTENSION_DESCRIPTION "DICOM plugin to import dynamic PET studies as Volume Sequences in 3D Slicer.")
set(EXTENSION_ICONURL "")
set(EXTENSION_SCREENSHOTURLS "")
set(EXTENSION_DEPENDS "NA")
set(EXTENSION_STATUS "Beta")

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
