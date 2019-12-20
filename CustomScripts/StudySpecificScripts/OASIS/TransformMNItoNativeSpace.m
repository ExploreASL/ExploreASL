%% Transform MNI to native space

MNInii = {fullfile(x.D.MapsDir,'Atlases','VascularTerritories','LabelingTerritories.nii')};

for iNii = 1:length(MNInii)
    [~, Ffile] = xASL_fileparts(MNInii{iNii});
    for iSubject=1:x.nSubjects
        PathOut = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, 'ASL_1', [Ffile '.nii']);
        InverseSpace = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, 'ASL_1', 'CBF.nii');
        DeformationPath = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, 'ASL_1', 'y_ASL.nii');
        AffineTrans = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, 'ASL_1', 'mean_PWI_Clipped_sn.mat');
        xASL_spm_deformations(x, MNInii{iNii}, PathOut, 0, InverseSpace, AffineTrans, DeformationPath);
    end
end