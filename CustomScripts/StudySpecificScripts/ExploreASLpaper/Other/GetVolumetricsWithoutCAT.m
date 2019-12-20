%% Get Volumetrics post-hoc (i.e. after CAT12 volume .mat files were deleted)

for iS=1:x.nSubjects
    xASL_TrackProgress(iS,x.nSubjects);
    C1File = fullfile(x.D.ROOT, x.SUBJECTS{iS}, 'c1T1.nii');
    C2File = fullfile(x.D.ROOT, x.SUBJECTS{iS}, 'c2T1.nii');
    C3File = fullfile(x.D.ROOT, x.SUBJECTS{iS}, 'c3T1.nii');

    if xASL_exist(C1File,'file') && xASL_exist(C2File,'file') && xASL_exist(C3File,'file')
%         C1nii = xASL_io_ReadNifti(C1File);
%         VoxelSize= prod(C1nii.hdr.pixdim(2:4));

        C1IM = xASL_io_Nifti2Im(C1File);
        C2IM = xASL_io_Nifti2Im(C2File);
        C3IM = xASL_io_Nifti2Im(C3File);

        ICV = sum(C1IM(:)) + sum(C2IM(:)) + sum(C3IM(:));

        GMvol = sum(C1IM(:)); % *VoxelSize;
        GM_ICV(iS,1)   = GMvol/ICV;
    end
end

IndexN = find(GM_ICV==min(GM_ICV));
x.SUBJECTS{IndexN}
