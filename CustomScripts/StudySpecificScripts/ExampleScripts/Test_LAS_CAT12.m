%% Loop over strength local adaptive segmentation
% PathOri = '/Users/henk/ExploreASL/ASL/Test_Dent_MCI-0013/analysis/MCI-0013/T1.nii';
% [Fpath, Ffile, Fext] = xASL_fileparts(PathOri);
% 
% for iStrength=1:20
%     LASstr = 0.05*iStrength;
%     
%     
%     NewDir = [Fpath '_T1_LAS' num2str(LASstr)];
%     xASL_adm_CreateDir(NewDir);
%     NewPath = fullfile(NewDir, 'T1.nii');
%     xASL_Copy(PathOri, NewPath);
%     
%     matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr = LASstr;
%     matlabbatch{1}.spm.tools.cat.estwrite.data{1,1} = NewPath;
%     
%     spm_jobman('run',matlabbatch);
% end

x =ExploreASL_Master('',0);

Rdir = '/Users/henk/ExploreASL/ASL/Test_Dent_MCI-0013/analysis';

for iStrength=1:20
    LASstr = 0.05*iStrength;
    NewDir = fullfile(Rdir, ['MCI-0013_T1_LAS' num2str(LASstr)]);
    if ~exist(NewDir,'dir')
        continue;
    end
    
    PathIn = {fullfile(NewDir, 'T1.nii'), fullfile(NewDir, 'mri', 'p1T1.nii')};
    PathOut = {fullfile(NewDir, 'rT1.nii') fullfile(NewDir, 'rc1T1.nii')};
    
    xASL_spm_deformations([], PathIn, PathOut, [], [], [], fullfile(NewDir, 'mri', 'y_T1.nii'));
    
    ImIn = {fullfile(NewDir, 'rT1.nii'), fullfile(NewDir, 'rc1T1.nii')};
    % slice 49
    IntScale = [0.75 1];
    NamePrefix = [];
    ColorMap = {x.S.gray, x.S.red};
    x.S.TraSlices = 49;
    bClip = [];
    MaskIn = [];
    bWhite = 0;
    MaxWindow = [];
    bTransparancy = [];
    
    ImOut{iStrength} = xASL_im_CreateVisualFig(x, ImIn, [], IntScale, NamePrefix, ColorMap, bClip, MaskIn, bWhite, MaxWindow, bTransparancy);
end

clear ImConcat
for iStrength=1:5
    if ~isempty(ImOut{iStrength})
        ImConcat(:,:,:,iStrength) = ImOut{iStrength};
    end
end

clear ImTotal
for iDim=1:size(ImConcat,3)
    ImTotal(:,:,iDim) = xASL_im_TileImages(squeeze(ImConcat(:,:,iDim,:)),size(ImConcat,4));
end
    
figure(1);imshow(ImTotal, 'InitialMagnification', 500)