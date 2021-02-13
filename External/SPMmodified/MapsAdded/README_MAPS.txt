The maps in this folder all originate from SPM and/or CAT12, and thus distributed under GNU GPL.

brainmask.nii: ICBM_152_nonlinear_symm brainmask, resampled to 1.5 mm MNI
brainmask_supratentorial.nii: same, but without structures outside cerebrum
CentralWM_QC.nii: ICBM_152_nonlinear_symm pWM eroded to center & resampled to 1.5 mm MNI
DeepWM.nii: ICBM_152_nonlinear_symm pWM eroded to deep WM & resampled to 1.5 mm MNI
Identity_Deformation_y_T1.nii: SPM deformation field, modified for no displacement, but resample to 1.5 mm MNI
Identity_sn.mat: SPM affine +DCT uniform deformation parameters, set for zero displacement
LabelColors.mat: iteration of different colors used to create labels
LeftRight.nii: 1.5 mm MNI template with 1 left and 2 right
ParenchymNarrow.nii: ICBM_152_nonlinear_symm wholeBrain (GM+WM) mask that is eroded & resampled to 1.5 mm MNI
rbrainmask.nii: ICBM_152_linear brainmask resampled to 1.5 mm MNI
rc1T1.nii & rc2T1.nii, : ICBM pGM & pWM templates resampled to 1.5 mm MNI
rc1T1_ASL_res.nii & rc2T1_ASL_res.nii: same but smoothed to typical ASL resolution (4x4x4mm)
rc3T1.nii & rc3T1_ASL_res.nii: created from rc1T1 & rc2T1, see below
rgrey.nii: SPM OldSeg pGM (ICBM_152_lin) resampled to 1.5 mm MNI
rT1.nii: ICBM T1 template resampled to 1.5 mm MNI
TotalGM.nii, TotalWM.nii, WholeBrain.nii: ICBM wholebrain mask, resampled to 1.5 mm MNI

MNI structural atlas, available at fsl.fmrib.ox.ac.uk/fsl/fslwiki/Atlases, but created under the ICBM license.

9 anatomical structural regions, kindly provided by Jack Lancaster at the Research Imaging Center, UTHSCSA, Texas (originally from the McConnell Brain Imaging Centre, MNI).

References:

Collins et al. Automatic 3-D model-based neuroanatomical segmentation. Human Brain Mapping 3(3): 190-208. (1995)
Mazziotta et al. A probabilistic atlas and reference system for the human brain: International Consortium for Brain Mapping (ICBM). Phil. Trans. Royal Soc. B Biol. Sci. 356(1412):1293-1322 (2001)


Several of these maps stem from the ICBM_152_linear & the ICBM_152_nonlinear atlas, which is free for all use (see license below):

Copyright (C) 1993–2009 Louis Collins, McConnell Brain Imaging Centre, Montreal Neurological Institute, McGill University. Permission to use, copy, modify, and distribute this software and its documentation for any purpose and without fee is hereby granted, provided that the above copyright notice appear in all copies. The authors and McGill University make no representations about the suitability of this software for any purpose. It is provided “as is” without express or implied warranty. The authors are not responsible for any data loss, equipment damage, property loss, or injury to subjects or patients resulting from the use or misuse of this software package. 

GhostSignalRatio.nii: template ROIs for the Ghost-To-Signal Ratio,
created from the ICBM152 brainmask.nii as above, per
REF: https://mriqc.readthedocs.io/en/stable/iqms/bold.html
1 = Signal
2 = Non-ghost
3 = Ghost

rc3T1.nii & rc3T1_ASL_res.nii: created by subtracting pGM and pWM from pMask, using the code below. While downloading a ICBM CSF map is probably nicer, this way would also need to update the rc1T1 & rc2T1, hence this way is more backward compatible/reproducible.

PathGM = '/Users/henk/ExploreASL/ExploreASL/External/SPMmodified/MapsAdded/rc1T1.nii';
PathWM = '/Users/henk/ExploreASL/ExploreASL/External/SPMmodified/MapsAdded/rc2T1.nii';
PathCSF = '/Users/henk/ExploreASL/ExploreASL/External/SPMmodified/MapsAdded/rc3T1.nii';
PathMask = '/Users/henk/ExploreASL/ExploreASL/External/SPMmodified/MapsAdded/brainmask.nii';

pGM=xASL_io_Nifti2Im(PathGM);
pWM=xASL_io_Nifti2Im(PathWM);
pMask = xASL_io_Nifti2Im(PathMask);

pCSF = pMask - min(pGM + pWM, ones(size(pGM)));
pCSF(pCSF<0) = 0;
pCSF(pCSF>1) = 1;

xASL_io_SaveNifti(PathGM, PathCSF, pCSF, [], 0);

% Same for rc3T3_ASL_res
PathGM = '/Users/henk/ExploreASL/ExploreASL/External/SPMmodified/MapsAdded/rc1T1_ASL_res.nii';
PathWM = '/Users/henk/ExploreASL/ExploreASL/External/SPMmodified/MapsAdded/rc2T1_ASL_res.nii';
PathCSF = '/Users/henk/ExploreASL/ExploreASL/External/SPMmodified/MapsAdded/rc3T1_ASL_res.nii';
PathMask = '/Users/henk/ExploreASL/ExploreASL/External/SPMmodified/MapsAdded/brainmask.nii';

pGM=xASL_io_Nifti2Im(PathGM);
pWM=xASL_io_Nifti2Im(PathWM);
pMask = xASL_io_Nifti2Im(PathMask);

pCSF = pMask - min(pGM + pWM, ones(size(pGM)));
pCSF(pCSF<0) = 0;
pCSF(pCSF>1) = 1;

xASL_io_SaveNifti(PathGM, PathCSF, pCSF, [], 0);