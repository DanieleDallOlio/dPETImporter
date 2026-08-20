# dPETImporter

<p align="center">
  <img src="dPETImporter.png" alt="dPETImporter" width="180">
</p>


| **Author** | **Project** |
|:----------:|:-----------:|
| [**Daniele Dall'Olio**](https://github.com/DanieleDallOlio) | dPETImporter |

**dPETImporter** is a DICOM import plugin for **3D Slicer** that recognizes
**dynamic PET (dPET)** studies and loads them as time-resolved **Volume Sequences**.

The plugin is designed to provide a reliable import layer for dynamic PET
workflows by reconstructing temporal frames from DICOM metadata and making the
result immediately available to Slicer's Sequences infrastructure.

dPETImporter is also distributed as part of the
[SlicerDynamicPET](https://github.com/DanieleDallOlio/SlicerDynamicPET)
extension.

---

## Features

- Automatic detection of **dynamic PET series**
- Loads dynamic datasets as **Slicer Volume Sequences**
- Supports both:
  - **3D-per-frame** datasets
  - **2D slice stacks per frame**
- Extracts and stores **frame timing metadata**
- Optional **SUVbw conversion during loading** when the required DICOM metadata are available
- Integrates directly with the **3D Slicer DICOM module**
- Creates the corresponding **Sequence Browser** for visualization and navigation

---

## Installation

### As part of SlicerDynamicPET

The recommended installation method is through the
**SlicerDynamicPET** extension.

Once the extension is installed and Slicer is restarted, dPETImporter registers
itself automatically with the DICOM module.

### Standalone development

The repository can also be built independently as a scripted Slicer extension
for development and testing.

---

## Usage

1. Open the **DICOM** module in 3D Slicer.
2. Import a dynamic PET DICOM study into the Slicer DICOM database.
3. Select the detected dynamic PET series.
4. Choose the dPETImporter load option when available.
5. Load the study as a **Volume Sequence**.
6. Navigate through the reconstructed temporal frames using the generated
   **Sequence Browser**.

---

## Output

The importer creates:

- a **Sequence Node** containing the reconstructed dynamic PET frames;
- a **Sequence Browser Node** for visualization and temporal navigation.

Frame timing and PET-related information are stored as MRML attributes where
available, so that downstream modules can access the acquisition timing
information needed for kinetic analysis.

---

## SUV Conversion

If sufficient metadata are available, dPETImporter can convert image values to
**SUVbw** during loading.

The required information includes, where applicable:

- radiopharmaceutical administration/start time;
- injected activity;
- radionuclide half-life;
- patient weight;
- decay-correction information;
- corrected-image flags.

SUVbw conversion is only applied when the required metadata are available and
the series is appropriately corrected, including attenuation and decay
correction.

If these requirements are not satisfied, the original image values are retained.

---

## Integration with SlicerDynamicPET

dPETImporter provides the DICOM import component of
[SlicerDynamicPET](https://github.com/DanieleDallOlio/SlicerDynamicPET).

The imported dynamic PET Volume Sequence can then be used by the DynamicPET
module for:

- time-activity curve extraction;
- graphical kinetic analysis;
- compartment-model fitting;
- voxel-wise parametric imaging.

The importer remains implemented as an independent Python component so that its
DICOM loading logic can be maintained and tested separately from the C++
kinetic-modeling module.

---

## Repository structure

```text
dPETImporter/
├── CMakeLists.txt
├── README.md
└── dPETImporter/
    ├── CMakeLists.txt
    ├── dPETImporter.py
    ├── dPETImporterPlugin.py
    └── reader.py
```

---

## Status

dPETImporter is research software under active development.

It is intended for research use and should be validated for the specific
datasets and workflows in which it is used.

---

## License

dPETImporter is distributed under the **MIT License**.

See [LICENSE](LICENSE) for details.

---

## Author

**Daniele Dall'Olio**  
University of Bologna

- [GitHub](https://github.com/DanieleDallOlio)
- [University profile](https://www.unibo.it/sitoweb/daniele.dallolio)

<a href="https://github.com/UniboDIFABiophysics">
<div class="image">
<img src="https://cdn.rawgit.com/physycom/templates/697b327d/logo_unibo.png" width="45" height="45">
</div>
</a>

---

## Acknowledgments

dPETImporter is developed for use within the 3D Slicer ecosystem.

Thanks to the 3D Slicer community and contributors for providing the DICOM,
Sequences, MRML, and Python infrastructure on which this plugin relies.
