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





