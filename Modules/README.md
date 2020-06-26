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
x = ExploreASL_Initialize(DataParPath, ProcessData, iWorker, nWorkers);
% Run the structural module
[~, x] = xASL_Iteration(x,'xASL_module_Structural');
```

----
## 2. Module ASL

#### Function

```matlab
[result, x] = xASL_module_ASL(x)
```

#### Description
This ExploreASL module processes the ASL images, i.e. ASL4D, M0, etc (if present), on an individual (i.e. session-to-session, where session indicates BIDS "run") basis. Both 2D and 3D options are automatially chosen, as well as processing of time-series (if available), such as motion correction and outlier exclusion.

#### Workflow

1. **TopUpASL**: FSL TopUp geometric distortion correction (if M0 images with reversed phase-encoding are present)
2. **RealignASL**: If time-series are present, motion correction and outlier exclusion (ENABLE)
3. **RegisterASL**: Registration of ASL to T1w anatomical images (if lacking, to MNI images)
4. **ResliceASL**: Resample ASL images to standard space
5. **PreparePV**: Create partial volume images in ASL space with ASL resolution
6. **ProcessM0**: M0 image processing
7. **Quantification**: CBF quantification
8. **CreateAnalysisMask**: Create mask using FoV, vascular outliers & susceptibility atlas
9. **VisualQC_ASL**: Generate QC parameters & images
10. **WADQC**:  QC for WAD-QC DICOM server (OPTIONAL)

#### Recommended usage

```matlab
[~, x] = xASL_Iteration(x,'xASL_module_ASL');
```

----
## 3. Module Population

#### Function

```matlab
[result, x] = xASL_module_Population(x)
```

#### Description
This ExploreASL module processes all available images on the group level. It assumes that all images were adequately processed in the previous modules.

#### Workflow

1. **CreatePopulationTemplates**: Create population average images, to compare scanners, cohorts etc without physiological variance
2. **CreateAnalysisMask**: Generate a group-level mask by combining individuals masks, for ROI-based analysis & VBA
3. **CreateBiasfield**: When there are multiple scanners, create scanner-specific biasfields (uses Site.mat for this)
4. **GetDICOMStatistics**: Create TSV file with overview of DICOM parameters
5. **GetVolumeStatistics**: Create TSV file with overview of volumetric parameters
6. **GetMotionStatistics**: Create TSV file with overview of motion parameters
7. **GetROIstatistics**: Create TSV file with overview of regional values (e.g. qCBF, mean control, pGM etc)
8. **SortBySpatialCoV**: Sort ASL_Check QC images by their spatial CoV in quality bins
9. **DeleteAndZip**: Delete temporary files and gzip all NIfTIs

#### Recommended usage

```matlab
[~, x] = xASL_Iteration(x,'xASL_module_Population');
```

