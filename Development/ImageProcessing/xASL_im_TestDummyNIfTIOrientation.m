% Example paths
FinalPath = 'ASLspaceButT1wOrientation.nii';
DummyPath = 'OrientationDummy.nii';
PathSrc = '../c1T1.nii';
PathRef = 'ASL4D.nii';

if xASL_exist(BackupPath)
    xASL_Copy(BackupPath, PathSrc, 1);
end
xASL_Copy(PathSrc, BackupPath, 1);
xASL_Copy(PathSrc, CopyPath,1);

% for example, let's apply only the dim order & zooming
bApplyZoom = 1;
bApplyRotation90 = 1;
bApplyRotationMinor = 0;
bApplyTranslation = 0;

xASL_im_DummyOrientationNIfTI(PathSrc, PathRef, DummyPath, bApplyRotationMinor, bApplyRotation90, bApplyZoom, bApplyTranslation);

%% Resample the source image into the new space/rotation
xASL_spm_reslice(DummyPath, PathSrc, [], [], 1, FinalPath);
xASL_delete(DummyPath);