c1File  = 'C:\Backup\ASL\Carolina\Trial_Lesion_nii\dartel\rc1T1_001_Lizama.nii';
c2File  = 'C:\Backup\ASL\Carolina\Trial_Lesion_nii\dartel\rc2T1_001_Lizama.nii';

c1FileNaN  = 'C:\Backup\ASL\Carolina\Trial_Lesion_nii\dartel\rc1T1_NaN_001_Lizama.nii';
c2FileNaN  = 'C:\Backup\ASL\Carolina\Trial_Lesion_nii\dartel\rc2T1_NaN_001_Lizama.nii';


LesionFile  = 'C:\Backup\ASL\Carolina\Trial_Lesion_nii\dartel\rLesion_T1_1_001_Lizama.nii';

c1IM    = xASL_io_Nifti2Im(c1File);
c2IM    = xASL_io_Nifti2Im(c2File);
LesionIM= xASL_io_Nifti2Im(LesionFile);

c1IM(logical(LesionIM))  = NaN;
c2IM(logical(LesionIM))  = NaN;

xASL_io_SaveNifti(c1File,c1FileNaN,c1IM);
xASL_io_SaveNifti(c2File,c2FileNaN,c2IM);
