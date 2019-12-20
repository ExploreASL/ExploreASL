x=ExploreASL_Master('',0);

DirP{1}  = 'C:\BackupWork\ASL\EPAD_xASL_Paper\analysis_Sleep\101'; % no atrophy, medium age
DirP{2}  = 'C:\BackupWork\ASL\EPAD_xASL_Paper\analysis_EPAD\040EPAD00040'; % slight atrophy fitting with age
DirP{3}  = 'C:\BackupWork\ASL\EPAD_xASL_Paper\analysis_EPAD\040EPAD00025'; % near-severest
DirP{4}  = 'C:\BackupWork\ASL\EPAD_xASL_Paper\analysis_Harmy\HD004_2'; % severest atrophy

for iL=1:length(DirP)
    %% Path definition - admin
    [RootDir{iL}, SubjName{iL}] = fileparts(DirP{iL});
    Path_PV_MNI{iL} = fullfile(RootDir{iL},'Population',['PV_pGM_' SubjName{iL} '.nii']);
    Path_PV{iL} = fullfile(DirP{iL},'PV_pGM.nii');

    Path_c1T1{iL} = fullfile(DirP{iL},'c1T1.nii');
    Path_c2T1{iL} = fullfile(DirP{iL},'c2T1.nii');
    Path_c3T1{iL} = fullfile(DirP{iL},'c3T1.nii');
    Path_T1{iL} = fullfile(DirP{iL},'T1.nii');
    Path_T1w{iL} = fullfile(DirP{iL},'T1_masked.nii');
    Path_T1w_MNI{iL} = fullfile(DirP{iL},'T1_masked_MNI.nii');
    Path_yT1{iL} = fullfile(DirP{iL},'y_T1.nii');
    Path_yASL{iL} = fullfile(DirP{iL},'ASL_1','y_ASL.nii');

    %% Create masked T1
    MaskIM{iL} = (xASL_io_Nifti2Im(Path_c1T1{iL})+xASL_io_Nifti2Im(Path_c2T1{iL})+xASL_io_Nifti2Im(Path_c3T1{iL}))>0.5;
    T1IM{iL} = xASL_io_Nifti2Im(Path_T1{iL});
    MaskIM{iL} = single(MaskIM{iL}).*T1IM{iL};
    xASL_io_SaveNifti(Path_c1T1{iL}, Path_T1w{iL}, MaskIM{iL}, [], 0);

    %% Reslice/re-center to MNI
    xASL_spm_reslice(x.D.ResliceRef, Path_T1w{iL}, x, Path_T1w_MNI{iL}, 1);

    %% Warp PV image back to T1w resolution
    xASL_spm_deformations(x, Path_PV_MNI{iL}, Path_PV{iL}, 1, Path_T1w_MNI{iL}, [], Path_yASL{iL});

    %% Open PV & threshold
    PVim{iL} = xASL_io_Nifti2Im(Path_PV{iL});
    Im50{iL} = PVim{iL}>0.5.*max(PVim{iL}(:));
    Im70{iL} = PVim{iL}>0.7.*max(PVim{iL}(:));
end

%% Construct Figure
TraSlices = [56 62 57 59];
bWhite = true;
clear NativeT1w Im50im Im70im RowIM T1im MaskIn
for iL=1:length(DirP)
    x.S.TraSlices = TraSlices(iL);
    T1im{iL} = Path_T1w_MNI{iL};
    MaskIn{iL} = T1im{iL}>10;
    NativeT1w{iL} = xASL_im_CreateVisualFig(x, T1im{iL}, [], 1, [], x.S.gray, [], MaskIn{iL}, bWhite);
    Im50im{iL} = xASL_im_CreateVisualFig(x, {T1im{iL} Im50{iL}}, [], [1 0.65], [], {x.S.gray x.S.red}, [], {MaskIn{iL} MaskIn{iL}}, bWhite);
    Im70im{iL} = xASL_im_CreateVisualFig(x, {T1im{iL} Im70{iL}}, [], [1 0.65], [], {x.S.gray x.S.red}, [], {MaskIn{iL} MaskIn{iL}}, bWhite);
    RowIM{iL} = [NativeT1w{iL} Im50im{iL} Im70im{iL}];
end
%
figure(1);imshow([RowIM{1};RowIM{4}],'InitialMagnification',250,'border', 'tight')
