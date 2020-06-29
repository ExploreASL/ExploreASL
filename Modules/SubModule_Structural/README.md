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

This submodule aims to clean up WMH under- or oversegmentations in a conservatively & robust way, i.e. erring on the side of caution.
It uses input from the tissue class segmentation (e.g. CAT12) to repair the WMH segmentation (e.g. LST LPA/LGA or externally provided).
Note that before running the tissue segmentation, the T1w was (conservatively) filled for WMH lesions.
This function is not tested a lot, so mainly conservatively set up to improve the WMH volumetrics, rather than improve the registration.

This submodule contains the following steps:
0. Administration
1. Correct pGM islands inside pWM: WMH can have an intensity similar to GM on the T1w, which erroneously classifies them as GM instead of WM(H). The rule used here, is to define GM islands within the WM as clusters of pGM>0.05 for which 3 layers (dilations) have at least 95% pWM. For these islands, pGM is given 100% to pWM. 50% of pWM is given to pWMH (the pWMH/pNAWM distinction is made later in the pipeline, here still pWM=pWMH+pNAWM). The reason is that not all low T1w intensities within the WM are WMH, we still expect some lacunes, perivascular (Virchow-Robin) spaces, which could be considered pNAWM rather than pWMH.
2. Perform brainmasking & join masks
3. Correct any WMH inside GM or CSF -> here we assume that CAT12 did a good segmentation job. If pGM is larger than pWM & larger than pWMH, we consider a voxel to be pGM and remove the pWMH. This effectively removes pWMH segmentation noise in the GM or CSF, it doesn't correct any significant misclassification of WMH in the GM or CSF. If the WMH segmentation does a significant misclassification (e.g. setting the pWMH inside GM or CSF to a probability higher than GM or CSF is by tissue segmentation), this is lesion filled after the WMH segmentation, on the T1w, hence the tissue segmentation won't have a chance to correct this. Fortunately, most oversegmentations in the GM/CSF have low pWMH, as WMH segmentation algorithms already perform a light tissue prior-based clean up themselves.
4. Saving & file management
5. Prepare visuals for visual QC & file management

### Recommended usage

### Interface definition







