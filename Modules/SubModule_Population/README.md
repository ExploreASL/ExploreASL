# Submodules of the Population Module

----
## Delete Many Temp Files

### Function

```matlab
function xASL_adm_DeleteManyTempFiles(x)
```

### Description

This function removes as many files as possible.

### Workflow

### Recommended usage

### Interface definition


----
## Gzip All Files

### Function

```matlab
function xASL_adm_GzipAllFiles(ROOT, bFolder, bUseLinux)
```

### Description

This function zips all .nii files recursively into nii.gz. Note that the linux option cannot deal with spaces in filenames or directories.

### Workflow

### Recommended usage

### Interface definition



----
## Create Analysis Mask

### Function

```matlab
function [x] = xASL_im_CreateAnalysisMask(x, Threshold)
```

### Description

Creates population analysis mask.

### Workflow

This function takes the mean population-based probability maps of masks, thresholds and combines them:
1. Creation GM, WM & WholeBrain masks by p>0.5
2. Create, combine & save vascular, susceptibity & FoV masks:
    * MaskVascular
    * MaskSusceptibility = MaskSusceptibility & MaskFoV
3. Create & save VBA mask
    * MaskAnalysis = MaskVascular & MaskSusceptibility
    * x.S.VBAmask = MaskAnalysis & GMmask
4. Visualization: Creates a figure with columns being following maps/masks overlaid over mean population T1w:
    * FoV probability 0-50% missing voxels
    * Vascular 0-7.5% missing voxels
    * Susceptibility 0-50% missing voxels
    * Analysis mask

### Recommended usage

### Interface definition



----
## Sort By Spatial CoV

### Function

```matlab
function xASL_qc_SortBySpatialCoV(x, Threshold1, Threshold2)
```

### Description

Compute statistics for each ROI.

### Workflow

This function organizes the ASL QC images in **//analysis/Population/ASLCheck** into CBF, vascular and artifactual contrast per the spatial CoV thresholds defined above, in folders:

* //analysis/Population/ASLCheck/1_CBFContrast
* //analysis/Population/ASLCheck/2_VascularContrast
* //analysis/Population/ASLCheck/3_ArtifactualContrast

Invalid spatial CoV numbers (e.g. NaN) will go to:

* //analysis/Population/ASLCheck/4_Unknown_sCoV

Note: these outputs need to be visually checked; but pre-sorting them by spatial CoV puts them already in categories that are quick to skim through and manually correct by moving the images

The idea is then that only category 1 images are used for perfusion (CBF) analyses, both categories 1 & 2 for the vascular (sCoV) analyses, and the category 3 images excluded from analysis.

**PM**: this code does not include multiple sessions per subject yet!

**NB**: this code uses the **//analysis/Population/Stats/CoV_qCBF*TotalGM*.csv** file, make sure that this file isn't edited! 

### Recommended usage

### Interface definition



----
## Compute Ws CV

### Function

```matlab
function x = xASL_stat_ComputeWsCV(x)
```

### Description

Calculates the within and between-subject coefficient of variance (wsCV and bsCV respectively), to estimate the power to detect effects. This requires 4D images that have been split.

### Workflow

### Recommended usage

### Interface definition



----
## Get Acquisition Time

### Function

```matlab
function x = xASL_stat_GetAcquisitionTime(x)
```

### Description

Summarize the acquisition time of scans.

### Workflow

This functions collects the DICOM field AcquisitionTime from each json sidecar (& parms.mat for backward compatibility) and saves them in the participants.tsv. Additionally, it creates a AcquisitionTime histogram of the full study, which can be useful to check time of scanning -> can influence physiological CBF variability.

1. Collect times
2. Save times
3. Create time histogram

### Recommended usage

### Interface definition



----
## Get DICOM Statistics

### Function

```matlab
function xASL_stat_GetDICOMStatistics(x, ScanType, HasSessions, bOverwrite)
```

### Description

Collect DICOM metadata in CSV file to check.

### Workflow

This functions prints DICOM metadata (e.g. parameters used for quantification) and collects them in a single tsv (per BIDS).
Summarizes this for the total population. Can be useful to detect software upgrades, where only slight parameter changes can hint on quantification changes.

This function carries out the following steps:

1. Load & save individual parameter files
2. Print summary
3. Write TSV file

### Recommended usage

### Interface definition



----
## Get Motion Statistics

### Function

```matlab
function xASL_stat_GetMotionStatistics(x)
```

### Description

Summarize motion values.

### Workflow

This functions collects motion stats, with the following steps:

1. Collect motion data
2. If no data, skip this function
3. Print motion vs exclusion overview
4. Add motion data to participants.tsv

### Recommended usage

### Interface definition



----
## Get Volume Statistics

### Function

```matlab
function xASL_stat_GetVolumeStatistics(x)
```

### Description

Summarize volume values.

### Workflow

This functions collects motion stats, with the following steps:

1. Collect structural volume data
2. Collect WMH data
3. Add stats in participants.tsv

### Recommended usage

### Interface definition



----
## Create Biasfield

### Function

```matlab
function xASL_wrp_CreateBiasfield(x)
```

### Description

Create sequence-specific intensity biasfields.

### Workflow

This function creates a smooth biasfield as intensity map to normalize intensities between different sequences/scanner/sites within a single study.
This is a simple pragmatic approach and is not validated, but is the best we can do.

First acquires average additive & multiplicative factors for total GM, then does a smooth voxel-wise rescaling. 
This doesn't make an assumption whether site or sequence differences are additive or multiplicative, but rather fits them both. 
Global scaling it performed to GM CBF == 60 mL/100g/min

**NB**: make sure that sequence resolution differences have been taken in account before creating these biasfields.

**PM**: add normalization of between-subjects SD as well.

**PM**: are there other things we can normalize?

### Recommended usage

### Interface definition



----
## Create Population Templates

### Function

```matlab
function xASL_wrp_CreatePopulationTemplates(x, SaveUnmasked, bCompute4Sets, SpecificScantype, bSkipMissingScans, bRemoveOutliers, FunctionsAre)
```

### Description

ExploreASL Population module wrapper, creates population parametric images for each ScanType.

### Workflow

This function creates simple parametric images, a.k.a. templates, for different image/scan types, on population level, as well as for different sets (e.g. sites/scanners/cohorts, etc) if specified. By default these images are masked, and transformed into a single column, for quick computations with low memory usage. 
The default parametric images that are created are the mean, between-subject SD, and the maximal intensity projection (MIP). 
The latter can e.g. identify intra-vascular signal that is similar between different subjects. Other parametric maps can be decommented (now commented out for speed).

### Recommended usage

### Interface definition



----
## Get ROI statistics

### Function

```matlab
function xASL_wrp_GetROIstatistics(x)
```

### Description

Compute statistics for each ROI.

### Workflow

This wrapper organizes the computation of statistics for different ROIs in a \[1.5 1.5 1.5\] mm MNI space:

1. Load the atlas: xASL_stat_AtlasForStats
2. Organize TSV output name: using x.S.output_ID
3. Obtain the ROI statistics: xASL_stat_GetROIstatistics
4. Print statistics in TSV files: xASL_stat_PrintStats

### Recommended usage

### Interface definition



----
## Load 4D Mem Mapping Lesions ROIs

### Function

```matlab
function [x] = xASL_wrp_Load4DMemMapping_LesionsROIs(x)
```

### Description

Part of ExploreASL analysis module. Loads data & maps it to memory mapping file on disc, if not done before.

### Workflow

### Recommended usage

### Interface definition


