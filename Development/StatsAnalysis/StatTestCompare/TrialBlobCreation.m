%% Trial blob creation

temp        = xASL_io_ReadNifti('E:\Backup\ASL_E\KCL_stats\GE_trial\analysis\dartel\DARTEL_CBF_HART12_ASL_1.nii');

ASL         = temp.dat(:,:,80);

BLOB        = CreateBlob( 12, [size(ASL,1),size(ASL,2)], 1, xASL_stat_MeanNan(ASL(:))*4, 0 , [30 -10]);

STAT        = ASL+BLOB;

dip_image(xASL_im_rotate([ASL;BLOB;STAT],90))

Plot2Dsurface(BLOB, 'E:\Backup\ASL_E\KCL_stats\GE_trial\analysis\BlobTrial\BlobSurface.jpg')
Plot2Dsurface(ASL, 'E:\Backup\ASL_E\KCL_stats\GE_trial\analysis\BlobTrial\ASLSurface.jpg')
Plot2Dsurface(STAT, 'E:\Backup\ASL_E\KCL_stats\GE_trial\analysis\BlobTrial\ASLplusBLOBSurface.jpg')
