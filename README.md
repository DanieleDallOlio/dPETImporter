# dPETImporter

<p align="center">
  <img src="dPETImporter.png" alt="dPETImporter logo" width="180">
</p>

| **Author** | **Project** |
| :---: | :---: |
| [**Daniele Dall'Olio**](https://github.com/DanieleDallOlio) | dPETImporter |

**dPETImporter** is a DICOM import plugin for [3D Slicer](https://www.slicer.org/) that interprets supported **dynamic PET**, **static/whole-body PET**, and **DICOM Parametric Map** data. It preserves quantitative and acquisition metadata used by downstream PET workflows.

It is included as a Git submodule in [SlicerDynamicPET](https://github.com/DanieleDallOlio/SlicerDynamicPET) but is maintained as a separate Python-based component.

> **Status: Research software.** DICOM variants and reconstruction settings differ between scanners. Always validate frame timing, quantitative calibration, source-image geometry, and derived results before scientific interpretation. This software is not a medical device.

## Supported import paths

### Conventional dynamic PET

The dynamic PET loader supports recognized conventional PET DICOM organizations, including:

- **3D-per-frame:** one PET image volume corresponding to each temporal frame;
- **2D-slices-per-frame:** multiple DICOM image instances reconstructed as each temporal volume.

Where source metadata permit, the plugin identifies frame timing, reconstructs time-ordered Slicer **Volume Sequences**, and creates a corresponding **Sequence Browser** for display and navigation.

### Enhanced PET

The plugin also supports a validated subset of **Enhanced PET** multiframe DICOM, including temporal acquisitions represented in a **single SOP Instance**. It reconstructs individual PET temporal volumes and preserves mappings to the original source-image frame numbers, image geometry, and quantitative transformations.

Enhanced PET is not treated as an unrestricted catch-all format. The reader validates per-frame functional groups, temporal and spatial organization, supported quantitative mappings/units, and source geometry. Examples of unsupported cases include incompatible temporal interleaving or multiple spatial stacks, missing mandatory timing/geometry, and ambiguous or nonlinear quantitative mappings. An unsupported Enhanced PET organization may require an importer update rather than a manual reordering of data.

### Static and whole-body PET (optional loadable)

An optional import path loads supported **static or whole-body PET** as ordinary Slicer scalar volumes while preserving kinetic/acquisition metadata. It uses Slicer's scalar-volume loader for the image and adds PET timing, provenance, and quantitative information.

**Important:** The dPETImporter static PET DICOM loadable is **unchecked by default**. Select the loadable marked with `dPET kinetic metadata` when this metadata-preserving behavior is required, for example for the late-static component of a SlicerDynamicPET Multi-timepoint study.

The importer detects whether individual images have **uniform acquisition timing** or **spatially varying timing** (as can occur across whole-body bed positions). The timing distinction is preserved rather than silently assigning one global acquisition time. Spatially varying timing is not currently accepted as a single late-static observation by SlicerDynamicPET's Multi-timepoint TAC workflow.

### DICOM Parametric Map import

A separate dPETImporter DICOM plugin recognizes supported **DICOM Parametric Map (PM)** objects and loads their quantitative voxel data into Slicer volumes. The implementation uses the Python `highdicom` library for reading DICOM PM data and geometry.

The code can attempt to install `highdicom` when absent. Internet access or a compatible already-installed package may be needed. As with PET images, derived map values and units should be checked against the originating acquisition and analysis.

## Acquisition and quantitative metadata

Where available, the plugin persists information such as:

- acquisition type and dynamic-frame or static-acquisition timing;
- per-frame start/end times and durations;
- radionuclide, tracer, injection time, administered dose, and half-life;
- PET units, decay-correction convention, and SUVbw-related factors;
- source SOP Instance UIDs, Series/Study identities, frame-of-reference identity, and source-frame references;
- relevant metadata for the downstream Multi-timepoint and RTSTRUCT workflows.

Metadata are stored on relevant MRML nodes (for example, dynamic PET sequences or static PET volumes), including attributes such as `dPET.KineticMetadata` and `dPET.DICOM.FrameReferences` where applicable. Persisted provenance can support downstream export without fabricating source-image identities.

### SUVbw conversion

When enabled and scientifically valid, dPETImporter can convert supported activity-concentration values to **body-weight-normalized SUV (SUVbw)** during loading.

The conversion requires compatible PET units and sufficient radiopharmaceutical, patient, and decay-correction metadata. The importer distinguishes values already stored as SUVbw from Bq/mL data that require conversion, and retains native values when the conversion is not supported.

For static whole-body images with spatially varying timing, the importer deliberately avoids applying a potentially incorrect single whole-volume SUV factor. Check the `dPET.ValueType` and other quantitative metadata before assuming two acquisitions can be compared directly.

## Installation

### Within SlicerDynamicPET

The recommended integrated setup is a compatible build of [SlicerDynamicPET](https://github.com/DanieleDallOlio/SlicerDynamicPET), which includes a pinned dPETImporter submodule. When building SlicerDynamicPET from source, initialize its submodules:

```bash
git submodule update --init --recursive
```

Once the scripted module is installed/loaded in Slicer, its DICOM plugins register with Slicer's DICOM infrastructure.

### Standalone development

The repository may also be built as a scripted Slicer extension against a compatible Slicer development environment. Example:

```bash
git clone https://github.com/DanieleDallOlio/dPETImporter.git
cmake -S dPETImporter -B dPETImporter-build -DSlicer_DIR=/path/to/Slicer-build
cmake --build dPETImporter-build --config Release
```

A source checkout alone does not automatically install the DICOM plugins into an unrelated Slicer instance.

## Usage

### Dynamic PET

1. Import the PET DICOM study into the **Slicer DICOM database**.
2. Select the appropriate dPETImporter dynamic PET loadable.
3. Load the acquisition as a **Volume Sequence**.
4. Use the generated **Sequence Browser** to inspect the temporal frames.
5. Verify the number and order of frames, start/end times, durations, image units, decay correction, and source-image references.

### Static PET for Multi-timepoint analysis

1. Import the static PET series through Slicer's DICOM module.
2. Explicitly check the optional loadable named with `dPET kinetic metadata`.
3. Load the scalar volume and inspect its timing mode and quantitative metadata.
4. In SlicerDynamicPET, select this PET along with the earlier dynamic acquisition in **Multi-timepoint** mode.

A late-static image with spatially varying slice/bed timing is not a substitute for a uniformly timed static observation.

### Parametric Maps

Import a supported DICOM PM object through Slicer's DICOM module and review the resulting scalar volume and its quantitative metadata. A DICOM PM represents a **derived parametric image**, not an additional dynamic PET acquisition.

## Output and integration

Depending on the selected import path, dPETImporter creates or updates:

- a `vtkMRMLSequenceNode` containing dynamic PET frames and a `vtkMRMLSequenceBrowserNode`;
- scalar volume nodes for static PET and imported DICOM Parametric Maps;
- MRML attributes carrying timing, DICOM provenance, units, decay correction, and calibration metadata.

These nodes can be used by SlicerDynamicPET for regional TAC extraction, compatible Multi-timepoint analysis, graphical/compartment modeling, parametric imaging, and provenance-aware exports.

## Limitations and validation

- Support for dynamic PET timing and quantitative metadata remains scanner- and reconstruction-dependent.
- Not every Enhanced PET multiframe organization or real-world-value mapping is supported.
- Whole-body images with differing slice/bed acquisition times cannot safely be modeled as one instantaneous static observation.
- Source DICOM should be retained for reproducibility and for workflows that explicitly require access to original instances.
- Metadata-preserving import is not a replacement for independently validating image geometry and quantitative calibration.

Please report reproducible incompatibilities, including the Slicer/importer commit, scanner and reconstruction software, anonymized relevant DICOM metadata, expected frame organization/timing, and error logs. **Do not attach identifiable clinical DICOM files to public issues.**

## Repository structure

```text
dPETImporter/
|-- CMakeLists.txt
|-- LICENSE
|-- README.md
|-- dPETImporter.png
`-- dPETImporter/
    |-- CMakeLists.txt
    |-- dPETImporter.py
    `-- dPETImporterPlugin.py
```

The importer logic, optional static loadable, Enhanced PET handling, and DICOM PM import code are implemented principally in `dPETImporterPlugin.py`.

## License

dPETImporter is distributed under the **MIT License**. See [LICENSE](LICENSE).

## Author

**Daniele Dall'Olio**  
University of Bologna

- [GitHub](https://github.com/DanieleDallOlio)
- [University profile](https://www.unibo.it/sitoweb/daniele.dallolio)

## Acknowledgments

Thanks to the 3D Slicer community for its DICOM, Sequences, MRML, volume, and Python infrastructure. dPETImporter is developed for use in the 3D Slicer ecosystem and integrated with SlicerDynamicPET.
