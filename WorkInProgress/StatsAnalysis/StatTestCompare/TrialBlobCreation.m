% Copyright 2015-2024 ExploreASL (Works In Progress code)
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

%% Trial blob creation

temp        = xASL_io_ReadNifti('E:\Backup\ASL_E\KCL_stats\GE_trial\analysis\dartel\DARTEL_CBF_HART12_ASL_1.nii');

ASL         = temp.dat(:,:,80);

BLOB        = CreateBlob( 12, [size(ASL,1),size(ASL,2)], 1, xASL_stat_MeanNan(ASL(:))*4, 0 , [30 -10]);

STAT        = ASL+BLOB;

dip_image(xASL_im_rotate([ASL;BLOB;STAT],90))

Plot2Dsurface(BLOB, 'E:\Backup\ASL_E\KCL_stats\GE_trial\analysis\BlobTrial\BlobSurface.jpg')
Plot2Dsurface(ASL, 'E:\Backup\ASL_E\KCL_stats\GE_trial\analysis\BlobTrial\ASLSurface.jpg')
Plot2Dsurface(STAT, 'E:\Backup\ASL_E\KCL_stats\GE_trial\analysis\BlobTrial\ASLplusBLOBSurface.jpg')
