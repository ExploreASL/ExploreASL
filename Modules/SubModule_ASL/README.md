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

0. Create FoV mask (native & MNI spaces)
1. Detect negative vascular signal (native & MNI spaces, within pGM>0.5)
2. Detect peak vascular signal (native & MNI spaces, within pGM==80% percentile on ASL image)
3. Brainmasking & FoV-masking (A) native & B) MNI spaces): Add WM vascular parts back to the mask (defined as pWM>0.8) & remove extracranial signal. In the WM, negative or peak signal is more expected from noise rather than from intra-vascular signal, not many big vessels exist in the WM.
4. Save vascular masks
5. Create susceptibility mask (standard space only). Here, we combine manually segmented susceptibility artifact regions in which a population-based susceptibility probability map is created. This map is combined (i.e. taking the product) with the mean control & PWI intensity distribution in these regions. This product is thresholded with the average of the 75th percentile & 15% of the intensity (for a bit more robustness against individual variability in sinus sizes).
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

Submodule of ExploreASL ASL Module, to prepare PV maps on ASL resolution.

### Workflow

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





