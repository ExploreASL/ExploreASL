# Modules

----
## 1. Module Structural

#### Function

```matlab
function [result, x] = xASL_module_Structural(x)
```

####Description
This first ExploreASL module processes the structural images, i.e. high-resolution **T1w** and **FLAIR** (if present), on an individual (i.e. subject-to-subject) basis. If a **FLAIR** is present, this is processed first to obtain a **WMH** mask to fill the hypointense lesions on the **T1w**, before segmenting the **T1w**. For the **T1w** segmentation this module uses **CAT12** by default but if this fails it falls back to **SPM** after trying to optimize the image contrast.

#### Workflow

1. **LinearReg\_T1w2MNI**: Ensure the alignment of subjects' anterior commissure (AC) with the AC in MNI & apply this to all images

2. **LinearReg\_FLAIR2T1w**: Align the FLAIR (if present) with T1w

3. **FLAIR\_BiasfieldCorrection**: Perform a biasfield correction (if not performed  by LST in following steps)

4. **LST\_Segment\_FLAIR_WMH**: Segment WMH lesions on FLAIR (if present)

5. **LST\_T1w\_LesionFilling_WMH**: Use WMH segmentation to fill lesions on T1w

6. **Segment\_T1w**: Tissue segmentation on T1w

7. **CleanUpWMH\_SEGM**: Extra WMH cleanup of some over- and under-segmentation

8. **Resample2StandardSpace**: Clone all images to standard space

9. **GetVolumetrics**: Obtain whole-brain volumes of GM, WM, CSF, WMH

10. **VisualQC**: Obtain QC parameters & save QC Figures

11. **DoWADQCDC**: QC for WAD-QC DICOM server (optional)

#### Recommended usage

```matlab
% Initialize dataset
DataParPath = "ExamplePath";
ProcessData = true;
iWorker = 1;
nWorkers = 1;
x = ExploreASL_Initialize(DataParPath, ProcessData, iWorker, nWorkers)
% Run the structural module
[~, x] = xASL_Iteration(x,'xASL_module_Structural');
```

----
## 2. Module ASL

...
