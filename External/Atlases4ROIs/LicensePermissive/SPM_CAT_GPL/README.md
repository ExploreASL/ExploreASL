The *masks* in this folder all originate from SPM and/or CAT12, and thus distributed under GNU GPL.
These are used as atlases/ROIs.
To differentiate from /External/SPMmodified/MapsAdded, which have the same license and origin but are *maps* used for image processing in ExploreASL.

DeepWM.nii: ICBM_152_nonlinear_symm pWM eroded to deep WM & resampled to 1.5 mm MNI
Supratentorial_GM_WM.nii: same as ICBM_152_nonlinear_symm brainmask, resampled to 1.5 mm MNI (/External/SPM_modified/MapsAdded/brainmask.nii), but without structures outside the cerebrum. For which we excluded the following Mindboggle regions: 1 cerebellum exterior, 2 cerebellum WM, 10 ventral dentate nucleus, 11 basal forebrain, 12-14 cerebellar vermis. (see /CustomScripts/Atlases/Supratentorial.m; private ExploreASL repository)
TotalGM.nii, TotalWM.nii, WholeBrain.nii: ICBM wholebrain mask, resampled to 1.5 mm MNI