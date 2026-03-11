# dPETImporter

<a href="https://github.com/UniboDIFABiophysics">
<div class="image">
<img src="https://cdn.rawgit.com/physycom/templates/697b327d/logo_unibo.png" width="90" height="90">
</div>
</a>

| **Authors**  | **Project** |
|:------------:|:-----------:|
| [**D. Dall'Olio**](https://github.com/DanieleDallOlio)  |  dPETImporter  |

**dPETImporter** is a DICOM plugin for **3D Slicer** that enables automatic loading of **dynamic PET (dPET)** studies as time-resolved image sequences.

The plugin detects dynamic PET acquisitions from DICOM metadata and loads them as **Volume Sequences**, making them immediately usable for visualization and quantitative analysis in Slicer.

---

## Features

- Automatic detection of **dynamic PET series**
- Loads dynamic datasets as **Slicer Volume Sequences**
- Supports both:
  - **3D-per-frame** datasets
  - **2D slice stacks per frame**
- Extracts and stores **frame timing metadata**
- Optional **SUVbw conversion during loading** when DICOM metadata are available
- Fully integrated in the **Slicer DICOM module**

---

## Installation

1. Install **3D Slicer**
2. Install the extension containing **dPETImporter**
3. Restart Slicer

The plugin will automatically register itself in the **DICOM module**.

---

## Usage

1. Open the **DICOM module**
2. Import a **dynamic PET DICOM dataset**
3. Select the detected dynamic PET series
4. Load it as a **Volume Sequence**

You can then navigate frames using the **Sequences browser**.

---

## Output

The importer creates:

- a **Sequence Node** containing the dynamic frames
- a **Sequence Browser Node** for visualization

Frame timing and SUV-related information are stored as MRML attributes.

---

## SUV Conversion

If sufficient metadata are available, images can be converted to **SUVbw** during loading.

Required DICOM metadata include:

- Radiopharmaceutical start time
- Injected dose
- Radionuclide half-life
- Patient weight
- Decay correction type
- Corrected image flags

SUV conversion is applied only if the series is **attenuation and decay corrected**.

---

## Author

* **Daniele Dall'Olio** [git](https://github.com/DanieleDallOlio), [unibo](https://www.unibo.it/sitoweb/daniele.dallolio)

## Acknowledgments

Thanks goes to all contributors of this project.
