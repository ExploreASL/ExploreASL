# Submodules of the Structural Module

----
## Clean Up WMH SEGM

### Function

```matlab
function xASL_wrp_CleanUpWMH_SEGM(x)
```

### Description

Submodule of ExploreASL Structural Module, that cleans up under- and over-segmentations of WMH SEGM.

### Workflow

This submodule aims to clean up **WMH** under- or oversegmentations in a conservatively & robust way, i.e. erring on the side of caution.
It uses input from the tissue class segmentation (e.g. **CAT12**) to repair the WMH segmentation (e.g. **LST LPA/LGA** or externally provided).
Note that before running the tissue segmentation, the **T1w** was (conservatively) filled for **WMH** lesions.
This function is not tested a lot, so mainly conservatively set up to improve the **WMH** volumetrics, rather than improve the registration.

This submodule contains the following steps:

0. Administration
1. Correct **pGM** islands inside **pWM**: **WMH** can have an intensity similar to **GM** on the **T1w**, which erroneously classifies them as **GM** instead of **WM(H)**. The rule used here, is to define **GM** islands within the WM as clusters of **pGM**>0.05 for which 3 layers (dilations) have at least 95% **pWM**. For these islands, **pGM** is given 100% to **pWM**. 50% of **pWM** is given to pWMH (the pWMH/pNAWM distinction is made later in the pipeline, here still **pWM=pWMH+pNAWM**). The reason is that not all low **T1w** intensities within the **WM** are **WMH**, we still expect some lacunes, perivascular (Virchow-Robin) spaces, which could be considered pNAWM rather than pWMH.
2. Perform brainmasking & join masks
3. Correct any **WMH** inside **GM** or **CSF** -> here we assume that **CAT12** did a good segmentation job. If **pGM** is larger than pWM & larger than pWMH, we consider a voxel to be **pGM** and remove the pWMH. This effectively removes pWMH segmentation noise in the **GM** or **CSF**, it doesn't correct any significant misclassification of **WMH** in the **GM** or **CSF**. If the **WMH** segmentation does a significant misclassification (e.g. setting the **pWMH** inside **GM** or **CSF** to a probability higher than **GM** or **CSF** is by tissue segmentation), this is lesion filled after the WMH segmentation, on the **T1w**, hence the tissue segmentation won't have a chance to correct this. Fortunately, most oversegmentations in the **GM/CSF** have low **pWMH**, as **WMH** segmentation algorithms already perform a light tissue prior-based clean up themselves.
4. Saving & file management
5. Prepare visuals for visual **QC** & file management

### Recommended usage

### Interface definition

----
## FLAIR Biasfield Correction

### Function

```matlab
function xASL_wrp_FLAIR_BiasfieldCorrection(x)
```

### Description

Submodule of ExploreASL Structural Module, that performs a biasfield correction on T1w & applies it on the FLAIR.

### Workflow

### Recommended usage

### Interface definition

----
## Get Volumetrics

### Function

```matlab
function xASL_wrp_GetVolumetrics(x)
```

### Description

Submodule of ExploreASL Structural Module, that obtains volumes from the tissue segmentations (& FLAIR WMH segmentations if they exist).

### Workflow

### Recommended usage

### Interface definition


----
## LST Segment FLAIR WMH

### Function

```matlab
function xASL_wrp_LST_Segment_FLAIR_WMH(x, rWMHPath, WMHsegmAlg)
```

### Description

Submodule of ExploreASL Structural Module, that performs a biasfield correction on T1w & applies it on the FLAIR.

### Workflow

### Recommended usage

### Interface definition



----
## LST T1w Lesion Filling WMH

### Function

```matlab
function xASL_wrp_LST_T1w_LesionFilling_WMH(x, rWMHPath)
```

### Description

Submodule of ExploreASL Structural Module, that performs lesion filling on T1w based on WMH segmented on the FLAIR.

### Workflow

### Recommended usage

### Interface definition



----
## Linear Reg FLAIR2T1w

### Function

```matlab
function xASL_wrp_LinearReg_FLAIR2T1w(x, bAutoACPC)
```

### Description

Submodule of ExploreASL Structural Module, that aligns FLAIR with T1w.

### Workflow

### Recommended usage

### Interface definition



----
## Linear Reg T1w2MNI

### Function

```matlab
function xASL_wrp_LinearReg_T1w2MNI(x, bAutoACPC)
```

### Description

Submodule of ExploreASL Structural Module, that aligns T1w with MNI.

### Workflow

### Recommended usage

### Interface definition



----
## Resample 2 Standard Space

### Function

```matlab
function xASL_wrp_Resample2StandardSpace(x)
```

### Description

Submodule of ExploreASL Structural Module, that resamples all structural images & derivatives.

### Workflow

### Recommended usage

### Interface definition



----
## Segment T1w

### Function

```matlab
function [x] = xASL_wrp_SegmentT1w(x, SegmentSPM12)
```

### Description

Submodule of ExploreASL Structural Module, that segments 3D T1 (or T2) scan.

### Workflow

### Recommended usage

### Interface definition



----
## Visual QC Structural

### Function

```matlab
function xASL_wrp_VisualQC_Structural(x)
```

### Description

Submodule of ExploreASL Structural Module, that performs several visualizations for .

### Workflow

### Recommended usage

### Interface definition




