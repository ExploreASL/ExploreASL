# Modules

----
## 1. Import Module

----
### ExploreASL\_Import

#### Format

```matlab
ExploreASL_Import(imPar[, bCopySingleDicoms, bUseDCMTK, bCheckPermissions, bRunDCM2NII, bClone2Source])
```

#### Description

Import batch T1, FLAIR, DWI, fMRI, M0, ASL data from dicom 2 NIfTI.
Uses dcm2niiX for the conversion, and additionally collects important DICOM header data
and puts them in \_parms.mat and .json sidecars to be used with the ExploreASL pipeline.
This function takes any folder input, but the folder input should be
specified in the imPar definition in the ExploreASL\_ImportConfig.m (later
to be converted to e.g. a CSV file). Follow the steps below, for study "MyStudy" located on "//MyDisk":

1) Make sure you have your DICOM data. Export them from XNAT, download them, or whatsoever
Create a root folder with study ID name, and put the DICOMs in any structure in the raw folder within the study ID root folder

----
## 2. Structural Module

----
### xASL\_module\_Structural

#### Format

```matlab
[result, x] = xASL_module_Structural(x)
```

#### Description
This first ExploreASL module processes the structural
images, i.e. high-resolution T1w and FLAIR (if present), on an individual (i.e. subject-to-subject) basis.
If a FLAIR is present, this is processed first to obtain a WMH mask to fill the hypointense lesions on the T1w,
before segmenting the T1w. For the T1w segmentation this module uses CAT12
by default but if this fails it falls back to SPM after trying to
optimize the image contrast. This module has the following steps/submodules/wrappers:

- 010\_LinearReg\_T1w2MNI         - Ensure the alignment of subjects' anterior commissure (AC) with the AC in MNI & apply this to all images
- 020\_LinearReg\_FLAIR2T1w       - Align the FLAIR (if present) with T1w
- 030\_FLAIR\_BiasfieldCorrection - Perform a biasfield correction (if not performed  by LST in following steps)
- 040\_LST\_Segment\_FLAIR\_WMH     - Segment WMH lesions on FLAIR (if present)
- 050\_LST\_T1w\_LesionFilling\_WMH - Use WMH segmentation to fill lesions on T1w
- 060\_Segment\_T1w               - Tissue segmentation on T1w
- 070\_CleanUpWMH\_SEGM           - Extra WMH cleanup of some over- and under-segmentation
- 080\_Resample2StandardSpace    - Clone all images to standard space
- 090\_GetVolumetrics            - Obtain whole-brain volumes of GM, WM, CSF, WMH
- 100\_VisualQC                  - Obtain QC parameters & save QC Figures
- 110\_DoWADQCDC                 - QC for WAD-QC DICOM server (OPTIONAL)


----
## 3. ASL Module

----
### xASL\_module\_ASL

#### Format

```matlab
[result, x] = xASL_module_ASL(x)
```

#### Description
This ExploreASL module processes the ASL
images, i.e. ASL4D, M0, etc (if present), on an individual (i.e. session-to-session, where session indicates BIDS "run") basis.
Both 2D and 3D options are automatially chosen, as well as processing of time-series (if available), such as motion correction and outlier
exclusion. This module has the following submodules/wrappers:

- 010\_TopUpASL          - FSL TopUp geometric distortion correction (if M0 images with reversed phase-encoding are present)
- 020\_RealignASL        - If time-series are present, motion correction and outlier exclusion (ENABLE)
- 030\_RegisterASL       - Registration of ASL to T1w anatomical images (if lacking, to MNI images)
- 040\_ResliceASL        - Resample ASL images to standard space
- 050\_PreparePV         - Create partial volume images in ASL space with ASL resolution
- 060\_ProcessM0         - M0 image processing
- 070\_Quantification    - CBF quantification
- 080\_CreateAnalysisMask- Create mask using FoV, vascular outliers & susceptibility atlas
- 090\_VisualQC\_ASL      - Generate QC parameters & images
- 100\_WADQC             - QC for WAD-QC DICOM server (OPTIONAL)


----
## 4. Population Module

----
### xASL\_module\_Population

#### Format

```matlab
[result, x] = xASL_module_Population(x)
```

#### Description
This ExploreASL module processes all available images on the
group level. It assumes that all images were adequately processed in the
previous modules. It will perform the following group-wise processing and
checks:

- 010\_CreatePopulationTemplates - Create population average images, to compare scanners, cohorts etc without physiological variance
- 020\_CreateAnalysisMask        - Generate a group-level mask by combining individuals masks, for ROI-based analysis & VBA
- 030\_CreateBiasfield           - When there are multiple scanners, create scanner-specific biasfields (uses Site.mat for this)
- 040\_GetDICOMStatistics        - Create TSV file with overview of DICOM parameters
- 050\_GetVolumeStatistics       - Create TSV file with overview of volumetric parameters
- 060\_GetMotionStatistics       - Create TSV file with overview of motion parameters
- 065\_GetRegistrationStatistics - Create TSV file with overview of the registration statistics
- 070\_GetROIstatistics          - Create TSV file with overview of regional values (e.g. qCBF, mean control, pGM etc)
- 080\_SortBySpatialCoV          - Sort ASL\_Check QC images by their spatial CoV in quality bins
- 090\_DeleteAndZip              - Delete temporary files and gzip all NIfTIs


