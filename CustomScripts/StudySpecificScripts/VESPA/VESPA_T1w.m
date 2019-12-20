ROOT = 'D:\Backup\ASL_E\VESPA_BaselineComparison';
for ii=1:22
    clear File1 File2 File3 IM
    File1   = fullfile(ROOT,'GE',sprintf('%03d',ii),'T1.nii');
    File2   = fullfile(ROOT,'CompleteT1s',['T1' sprintf('%03d',ii) '.nii']);
    File3   = fullfile(ROOT,'PH',sprintf('%03d',ii),'T1.nii');

    IM      = xASL_io_ReadNifti(File1);
    IM      = IM.dat(:,:,:,1);
    xASL_io_SaveNifti( File1, File3, IM );

    xASL_Move(File1,File2);
end


ROOT = 'D:\Backup\ASL_E\VESPA_BaselineComparison';
for ii=1:22
    clear File1 File2 File3 IM
    File1   = fullfile(ROOT,'CompleteT1s',['T1' sprintf('%03d',ii) '.nii']);
    File2   = fullfile(ROOT,'PH',sprintf('%03d',ii),'T1.nii');

    IM      = xASL_io_ReadNifti(File1);
    IM      = IM.dat(:,:,:,2);
    xASL_io_SaveNifti( File1, File2, IM );

end
