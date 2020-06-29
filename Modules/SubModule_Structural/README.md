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
1. Correct **pGM** islands inside **pWM**: **WMH** can have an intensity similar to **GM** on the **T1w**, which erroneously classifies them as **GM** instead of **WM(H)**. The rule used here, is to define **GM** islands within the **WM** as clusters of **pGM**>0.05 for which 3 layers (dilations) have at least 95% **pWM**. For these islands, **pGM** is given 100% to **pWM**. 50% of **pWM** is given to pWMH (the **pWMH/pNAWM** distinction is made later in the pipeline, here still **pWM=pWMH+pNAWM**). The reason is that not all low **T1w** intensities within the **WM** are **WMH**, we still expect some lacunes, perivascular (Virchow-Robin) spaces, which could be considered pNAWM rather than pWMH.
2. Perform brainmasking & join masks
3. Correct any **WMH** inside **GM** or **CSF** -> here we assume that **CAT12** did a good segmentation job. If **pGM** is larger than pWM & larger than **pWMH**, we consider a voxel to be **pGM** and remove the pWMH. This effectively removes pWMH segmentation noise in the **GM** or **CSF**, it doesn't correct any significant misclassification of **WMH** in the **GM** or **CSF**. If the **WMH** segmentation does a significant misclassification (e.g. setting the **pWMH** inside **GM** or **CSF** to a probability higher than **GM** or **CSF** is by tissue segmentation), this is lesion filled after the WMH segmentation, on the **T1w**, hence the tissue segmentation won't have a chance to correct this. Fortunately, most oversegmentations in the **GM/CSF** have low **pWMH**, as **WMH** segmentation algorithms already perform a light tissue prior-based clean up themselves.
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

Submodule of ExploreASL Structural Module, that performs a biasfield correction on **T1w** & applies it on the **FLAIR**.

### Workflow

This submodule performs a biasfield correction on **T1w** and applies it on **FLAIR**. This can be useful, when there are large lesions on the **FLAIR** that hamper capturing the biasfield nicely on the **FLAIR** itself. In such cases, the biasfield of the T1w might be easier to obtain and should be the same as the **FLAIR**, provided they are scanned in the same scan session (i.e.g same scanner, same coil).

**BE CAREFUL**: this submodule assumes that the biasfields of the **T1w** and **FLAIR** are comparable, which is not the case when one of the two (or both) are already biasfield corrected.

### Recommended usage

### Interface definition

----
## Get Volumetrics

### Function

```matlab
function xASL_wrp_GetVolumetrics(x)
```

### Description

Submodule of ExploreASL Structural Module, that obtains volumes from the tissue segmentations (& **FLAIR WMH** segmentations if they exist).

### Workflow

This submodule computes the total volume for each of the tissue classes & stores them in a **TSV** file (per **BIDS**).
This is computed from the native space segmentation derivatives (**GM, WM & CSF**), from which the **ICV** & relative volumes can be calculated. 
This is performed for **CAT12** or **SPM12** (whichever was used), and optionally for a **WMH_SEGM**.

### Recommended usage

### Interface definition


----
## LST Segment FLAIR WMH

### Function

```matlab
function xASL_wrp_LST_Segment_FLAIR_WMH(x, rWMHPath, WMHsegmAlg)
```

### Description

Submodule of ExploreASL Structural Module, that performs a biasfield correction on **T1w** & applies it on the **FLAIR**.

### Workflow

This submodule runs the **LST WMH** segmentation, either with **LGA** or **LPA**.
**LPA** is the default choice, it outperforms LGA a bit, depending on the image quality. 
These algorithms perform optimally with **3T** images, with good contrast. Generally, **LPA** oversegments whereas **LGA** undersegments.
The **LPA** oversegmentation is corrected in a later submodule.
If a **WMH_SEGM** already existed, the **LST** is run quickly as dummy only, to be replaced by the original **WMH_SEGM** image. This function has the following parts:

1. Reslice **FLAIR** (& **WMH_SEGM**, if exists) to T1w
2. Define parameters for segmentation
3. Run the segmentation
4. Replace by already existing **WMH_SEGM**
5. File management
6. Remove NaNs from segmentations & fix image edges

### Recommended usage

### Interface definition



----
## LST T1w Lesion Filling WMH

### Function

```matlab
function xASL_wrp_LST_T1w_LesionFilling_WMH(x, rWMHPath)
```

### Description

Submodule of ExploreASL Structural Module, that performs lesion filling on **T1w** based on **WMH** segmented on the **FLAIR**.

### Workflow

This submodule runs the **LST** **WMH**-based **T1w** lesion filling, which should improve the registration & segmentation of the **T1w** by e.g. **CAT12/SPM12**. The **WMH** can be either segmented in the previous submodule by **LST LGA/LPGA** or provided externally.
Before lesion filling, we clean up the WMH segmentation, to make the lesion filling a bit more conservative. 
Sometimes the **WMH** segmentation oversegments inside the GM (as there can be hyperintensities on the **FLAIR**) & we don't want to lesion-fill these on the **T1w** (which would turn their intensities in intensities similar to **WM**, leading to misclassifications by the **T1w** segmentation).
Note that this is submodule only performs the lesion filling, and the clean up is also performed for the lesion filling only.
A more thorough **WMH** clean up (for e.g. **WMH** volumetrics) is performed later in the Structural module, using also the results from the T1w segmentation.

Note when changing the lesion filling here, **LST** lesion filling expects a probability map, doesnt work nicely with binary mask
This function runs the following steps:

1. File management
2. Clean up the WMH segmentation used for lesion filling
3. Run lesion filling
4. Correction of too much/erronous lesion filling
5. File management

### Recommended usage

### Interface definition



----
## Linear Reg FLAIR 2 T1w

### Function

```matlab
function xASL_wrp_LinearReg_FLAIR2T1w(x, bAutoACPC)
```

### Description

Submodule of ExploreASL Structural Module, that aligns **FLAIR** with **T1w**.

### Workflow

This submodule registers **FLAIR** linearly to the **T1w**.
The same transformation is applied to all other related scans (**FLAIR**-segmented lesions, **WMH** specifically or other lesions).
This is required to enable the application of T1w derivatives (e.g. transformations to standard space, tissue segmentation) for **FLAIR** and vice versa (e.g. **WMH** lesion-filling).

### Recommended usage

### Interface definition



----
## Linear Reg T1w2MNI

### Function

```matlab
function xASL_wrp_LinearReg_T1w2MNI(x, bAutoACPC)
```

### Description

Submodule of ExploreASL Structural Module, that aligns **T1w** with **MNI**.

### Workflow

This submodule registers T1w linearly to the center of MNI space, a.k.a. **ACPC** alignment.
The same transformation is applied to all other related scans (**ASL4D, M0, FLAIR,** etc.).
This facilitates **MNI**-based algorithms (e.g. **SPM**-based segmentation), and allows for visual **QC** with all images roughly in the same space. This submodule first clips high values that can bias the registration algorithm, then performs a center of mass-based **ACPC** alignment, and then several iterations of **SPM** coregistration.
Assuming that this submodule is run at the start of ExploreASL, all NIfTI orientation matrices are restored before running the registration.

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

This submodule resamples all structural images & their derivatives to standard space. 
It uses the transformation fields that were obtained previously in the Structural module, concatenates all transformations into a single transformation (if not already done) & applies the transformation with a single interpolation (either trilinear for low quality or probability maps, or 2nd order B-spline). 
Finally, it computes the Jacobian determinants (i.e. the derivative of the transformation field) to obtain a map of the volumetric effects of the transformation. 
This Jacobian map is multiplied with the standard space resampled images, to restore their (local & global) total volume. 
The sum of volumes in native & standard space are compared as **QC**.
This submodule is not only part of the structural module, but can be repeated when the transformation map is edited, e.g. after longitudinal registration or after creation of a group-wise template.

### Recommended usage

### Interface definition



----
## Segment T1w

### Function

```matlab
function [x] = xASL_wrp_SegmentT1w(x, SegmentSPM12)
```

### Description

Submodule of ExploreASL Structural Module, that segments **3D T1** (or **T2**) scan.

### Workflow

This submodule segments high resolution structural/anatomical scans into **GM/WM/CSF**/soft tissue/bone/air tissue classes.
It will save **GM/WM/CSF** in native space, and the transformation from native to standard space.
This transformation includes Geodesic Shooting/DARTEL for **CAT12**.

This submodule contains the following steps:

1. Administration
2. Extra segmentation options by Jan Petr
3. Segmentation using **CAT12**
    * If **CAT12** fails, it will be repeated with higher contrast, higher strength affine preprocessing & less biasfield regularization
    * If **CAT12** fails twice, it will be skipped & **SPM12** will be run
4. Segmentation using SPM12
5. File management **CAT12**
6. File management lesions
7. Resample lesions to standard space
    * For the lesion masking. MORE EXPLANATION NEEDED BY JAN
8. Manage flowfields
    * Smooth combination non-linear flowfield outside the lesion & uniform flowfield within the lesion
9. File management

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

This submodule performs several visualizations for visual & quantitative **QC**.

1. After initial admin
2. It starts with the SPM UP parameters (courtesy of Cyril Pernet, his SPM UP scripts were made more robust & accurate by Jan & Henk, & are implemented here for **T1w** (& optionally **FLAIR**).
3. Then it performs a collection of visualizations
4. Also repeated specifically for lesions & manually provided **ROIs**
5. Finally, this contains a report of all missing raw & derivative files, in native & standard space, printing the NIfTI orientation matrix content before (**hdr.mat0**) & after registrations (**hdr.mat**). The determinant of these matrices should be the same, otherwise LeftRight has flipped. This should also be the same across a group scanned at the same scanner. Then various other **QC** functions are called & all are summarized in a PDF report.

### Recommended usage

### Interface definition




