# Submodules of the ASL Module

----
## Create Analysis Mask

### Function

```matlab
function xASL_wrp_CreateAnalysisMask(x)
```

### Description

Create analysis mask from a combination of FoV & removal of high and negative intravascular ASL voxels.

### Workflow

0. Create **FoV** mask (native & **MNI** spaces)
1. Detect negative vascular signal (native & **MNI** spaces, within **pGM**>0.5)
2. Detect peak vascular signal (native & MNI spaces, within **pGM**==80% percentile on **ASL** image)
3. Brainmasking & **FoV**-masking (A) native & B) **MNI** spaces): Add **WM** vascular parts back to the mask (defined as pWM>0.8) & remove extracranial signal. In the **WM**, negative or peak signal is more expected from noise rather than from intra-vascular signal, not many big vessels exist in the **WM**.
4. Save vascular masks
5. Create susceptibility mask (standard space only): Here, we combine manually segmented susceptibility artifact regions in which a population-based susceptibility probability map is created. This map is combined (i.e. taking the product) with the mean control & PWI intensity distribution in these regions. This product is thresholded with the average of the 75th percentile & 15% of the intensity (for a bit more robustness against individual variability in sinus sizes).
6. Create standard space CBF_masked image to visualize masking effect

### Recommended usage

### Interface definition


----
## Prepare PV

### Function

```matlab
function x = xASL_wrp_PreparePV(x, bStandardSpace)
```

### Description

Submodule of ExploreASL **ASL** Module, to prepare **PV** maps on **ASL** resolution.

### Workflow

This submodule prepares partial volume correction (**PVC**) by creating correct **PV** maps in **ASL** resolution, in native space, as well as in standard space if requested (to perform **PVC** in standard space).

**If** bStandardSpace:
1. Create dummy upsampled ASL scan, for registration
2. Reslice **pGM** & **pWM** to hi-res **ASL**
3. Estimate effective spatial resolution of **ASL**
4. Smooth **pGM** & **pWM** to this spatial resolution
5. Move smoothed tissue posteriors to **MNI** space

**Else:**
Run step 3 only, which will use the effective spatial resolution that is default for the respective sequence:
* **2D EPI**: \[1 1 1\] * VoxelSize
* **3D GRASE**: \[1.1 1.1 1.38\] * VoxelSize
* **3D spiral**: \[4.3 4.4 10.1\] * VoxelSize (assuming **GE** uses the upsampled 2x2x4 mm 
& run steps 1&2, but in native space these entail presmoothing & downsampling.)

### Recommended usage

### Interface definition


----
## Process M0

### Function

```matlab
function xASL_wrp_ProcessM0(x)
```

### Description

Submodule of ExploreASL ASL Module, for M0 image processing.

### Workflow

This submodule performs the image processing and quantification of **M0** maps (if they exist), with the following steps:

1. Register **M0** to mean control if it exists: Before registration, contrast is equalized between the images & biasfields are removed.
2. Quantify **M0** (correction scale slope & incomplete T1 recovery)
3. Masking & smoothing of **M0** image, either using:
* **A)** traditional technique (very sharp masking & little smoothing)
* **B)** new ExploreASL-specific technique:
* extrapolating outside mask (avoiding artifacts from too much or too little masking)
* smooth very extensively, to create a biasfield (increases robustness & comparison of M0 between sequences/patients)

Any **M0** will be processed here. Even if part of the subjects does not have an **M0**, since this can be later imputed, or an average population **M0** image could be used. Also, without background suppression and without an **M0**, the MeanControl image is before saved as **M0**, and will be processed here as well.

Note that any voxel-size differences between **M0** and ASL are allowed here: step 0B below rescales the **PD** inside an **M0** voxel to the same as the **ASL** resolution (assuming a voxel with half volume contains half the amount of protons). The M0 is further processed in standard space, and reduced to a biasfield. For the quantification in standard space, the **PWI** and **M0** are now by definition in the same space. Also, the standard space **M0** biasfield is resampled to the native **PWI** space (at the end of step **3B)**, ensuring that both are also in the same native space.

### Recommended usage

### Interface definition


----
## Quantify

### Function

```matlab
function xASL_wrp_Quantify(x, PWI_Path, OutputPath, M0Path, SliceGradientPath)
```

### Description

Submodule of ExploreASL ASL Module, that performs quantfication.

### Workflow

This submodule converts PWIs to quantified CBF maps (or related derivatives):

1. Load **PWI**
2. Prepare **M0**
3. Hematocrit & blood **T1** correction
4. **ASL** & **M0** parameters comparisons
5. Load SliceGradient
6. Initialize & define quantification parameters
7. Define labeling efficiency
8. Perform quantification
9. Save files
10. Perform FEAST quantification (if exist)

### Recommended usage

### Interface definition


----
## Realign ASL

### Function

```matlab
function xASL_wrp_RealignASL(x,bSubtraction)
```

### Description

Submodule of ExploreASL ASL Module, that realigns.

### Workflow

This submodule estimates motion by spm_realign, which uses a rigid-body registration (3 translations, 3 rotations). It runs ENABLE to reject outliers and provides a visualization. ENABLE, QC and visualizations are based on the Net Displacement Vector (NDV) (in mm): according to Pythagorean/Euclydian RMS.

[jiscmail](https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind1211&L=fsl&P=R34458&1=fsl&9=A&J=on&d=No+Match%3BMatch%3BMatches&z=4)

view this link for image of rotation roll, pitch and yaw [grcnasa](https://www.google.nl/search?q=rotation+pitch+yaw+roll&espv=2&tbm=isch&imgil=LW3Nn1K-L6Oc7M%253A%253B-aSyykkRityJoM%253Bhttp%25253A%25252F%25252Fwww.grc.nasa.gov%25252FWWW%25252Fk-12%25252Fairplane%25252Frotations.html&source=iu&usg=__MlLQ5VuyRbm6kZP0vBJlPxmfbkw%3D&sa=X&ei=TWfjU4WcK4bqyQPqu4Fo&ved=0CD8Q9QEwBQ&biw=1680&bih=946#facrc=_&imgdii=_&imgrc=LW3Nn1K-L6Oc7M%253A%3B-aSyykkRityJoM%3Bhttp%253A%252F%252Fwww.grc.nasa.gov%252FWWW%252Fk-12%252Fairplane%252FImages%252Frotations.gif%3Bhttp%253A%252F%252Fwww.grc.nasa.gov%252FWWW%252Fk-12%252Fairplane%252Frotations.html%3B709%3B533)

This submodule performs the following steps:
1. Estimate motion
2. Calculate and plot position and motion parameters
3. Threshold-free spike definition (based on ENABLE, but with t-stats rather than the threshold p<0.05)
4. Remove spike frames from nifti

### Recommended usage

### Interface definition


----
## Register ASL 

### Function

```matlab
function xASL_wrp_RegisterASL(x)
```

### Description

Submodule of ExploreASL ASL Module, that registers ASL to T1w (or potentially other structural images).

### Workflow

### Recommended usage

### Interface definition


----
## Resample ASL

### Function

```matlab
function xASL_wrp_ResampleASL(x)
```

### Description

Submodule of ExploreASL ASL Module, that reslices native space images to standard space.

### Workflow

### Recommended usage

### Interface definition


----
## Visual QC ASL

### Function

```matlab
function xASL_wrp_VisualQC_ASL(x)
```

### Description

Submodule of ExploreASL ASL Module, that performs several visualizations for QC.

### Workflow

### Recommended usage

### Interface definition





