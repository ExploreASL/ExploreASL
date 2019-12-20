These new functions were added to this modified SPM version, for backward compatibility,
improved reproducibility, nii.gz file handling, cost function masking, and allowing a low and high quality setting for quick testing. Graphical progress tracking replaced by non-graphical progress tracking (CLI). All under GNU GPL with SPM.



FILE HANDLING:
xASL_adm_ConvertSlash.m - manage different slash between linux/windows
xASL_adm_CreateDir.m - create new folder, including parent folders, if not already existing
xASL_adm_GetFileList - wrapper around SPM_select
xASL_adm_ManageMoCoMat - manage the external .mat file for timeseries
xASL_Copy - wrapper around Matlabs copy
xASL_Move - same as previous but moving instead of copying
xASL_delete - wrapper around Matlabs delete (manage .nii.gz, delete only if exists)
xASL_exist - wrapper around Matlabs exist, manage .nii.gz
xASL_fileparts - wrapper around Matlabs fileparts, manage .nii.gz
xASL_io_ReadNifti - wrapper around SPMs nifti, with .nii.gz support and extra checks
xASL_io_Nifti2im - Same as previous but immediately convert to image matrix
xASL_io_SaveNifti - store image matrix to NIfTI file
xASL_round - wrapper around Matlabs round, to allow determining nr floating points
xASL_spm_admin - manage SPM input
xASL_spm_deformations - wrapper around SPM deformations, specific combinations for ExploreASL 							toolbox
xASL_spm_reslice - wrapper around SPM reslice, manage affine, allowing different quality, etc
xASL_spm_smooth - wrapper around SPM smooth
xASL_SysCopy - part of xASL_Copy, without .nii.gz support (use when dealing with only .nii OR 					.nii.gz)
xASL_SysMove - same as previous but moving instead of copying
xASL_TrackProgress - print a percentage tracker on the screen, to follow progress without GUI
xASL_wrp_CreateASLDeformationField - downsample/smooth T1w flowfield to ASL resolution



IMAGE PROCESSING:
xASL_im_ConvertMap2Mask - in the name
xASL_im_ExtrapolateOverNaNs - fill all NaNs within image by smoothing with extrapolation
xASL_FixEdgesFlowField - registration/interpolating can give artifacts at flowfield edges, fix 							these
xASL_im_ndnanfilter - handle NaNs when smoothing
xASL_im_ResampleIM - resample images using Matlabs interp function
xASL_init_DefaultEffectiveResolution - manage ASL effective resolution



COST FUNCTION MASKING
xASL_im_LesionRemoval4CAT
xASL_im_SaveOriginal4CAT - save original flow field
xASL_wrpDARTELSaveIntermedTrans - save intermediate transformations while running DARTEL
xASL_wrp_GSSaveIntermedTrans - same as previous but for Geodesic Shooting (GS) instead of DARTEL



.NII.GZ HANDLING:
xASL_adm_UnzipNifti
xASL_adm_ZipFileNameHandling



INCREASE SPM REPRODUCIBILITY OVER MATLAB VERSIONS BY NEW CONVOLUTION
xASL_im_conv3Dsep(.m|c|mex|mex64|mexmaci64|mexw64)
xASL_mex_compile_all -> code to quickly compile these





Thomas Nichols, adapted from "Johns Gems" - CORRCLUSTH
The CorrClusTh.txt script is from John's Gems/Thom Nichols, and falls under the GNU GPL as well.