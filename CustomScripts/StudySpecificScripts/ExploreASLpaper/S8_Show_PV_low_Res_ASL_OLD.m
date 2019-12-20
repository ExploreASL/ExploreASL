x=ExploreASL_Master('',0);

DirP{1}  = 'C:\BackupWork\ASL\EPAD_xASL_Paper\analysis_EPAD\040EPAD00040'; % no atrophy
C:\BackupWork\ASL\EPAD_xASL_Paper\analysis_EPAD\040EPAD00025 % near-severest
DirP{2}  = 'C:\BackupWork\ASL\Harmy\analysis_Harmy\HD004_2'; % severest atrophy
DirP{3}  = 'C:\Users\henkj\Google Drive\WorkGeneral\PapersPending\01_Explore_ASL_paper\Figures\SupplFig8\pvGM_left_in_atrophy\Atrophy_vs_pGM_voxels_3Dspiral_pseudo_2DEPI_PVmaps\dartel';


for iL=1:2
    Path_T1w{iL}        = fullfile(DirP{iL},'T1_masked.nii');
    Path_T1w_MNI{iL}    = fullfile(DirP{iL},'T1_masked_MNI.nii');
    % Reslice/re-center to MNI
    xASL_spm_reslice( x.D.ResliceRef, Path_T1w{iL}, [], [], x.Quality, Path_T1w_MNI{iL}, 1 );

    IM_T1w_MNI{iL}      = xASL_io_Nifti2Im(Path_T1w_MNI{iL});
end



%% Created thresholded pGM images in ASL resolution, move to native space

Path_PV{1}         = fullfile(DirP{3}, 'PV_pGM_NoAtrophy.nii');
Path_PV{2}         = fullfile(DirP{3}, 'PV_pGM_Atrophy.nii');

for iL=1:2
    Path_PV_50_MNI{iL}  = fullfile(DirP{iL}, 'PV_pGM_50_MNI.nii');
    Path_PV_70_MNI{iL}  = fullfile(DirP{iL}, 'PV_pGM_70_MNI.nii');
    Path_PV_50{iL}      = fullfile(DirP{iL}, 'PV_pGM_50.nii');
    Path_PV_70{iL}      = fullfile(DirP{iL}, 'PV_pGM_70.nii');
    yPath{iL}           = fullfile(DirP{iL}, 'y_T1.nii');

    xASL_io_SaveNifti( Path_PV{iL}, Path_PV_50_MNI{iL}, single(xASL_io_Nifti2Im(Path_PV{iL})>0.5), 32, 0);
    xASL_io_SaveNifti( Path_PV{iL}, Path_PV_70_MNI{iL}, single(xASL_io_Nifti2Im(Path_PV{iL})>0.7), 32, 0);
    xASL_spm_deformations(x,Path_PV_50_MNI{iL},Path_PV_50{iL},0,Path_T1w_MNI{iL},[],yPath{iL});
    xASL_spm_deformations(x,Path_PV_70_MNI{iL},Path_PV_70{iL},0,Path_T1w_MNI{iL},[],yPath{iL});
%     xASL_spm_reslice( x.D.ResliceRef, Path_PV_50{iL}, [], [], x.Quality, Path_PV_50{iL}, 1 );
%     xASL_spm_reslice( x.D.ResliceRef, Path_PV_70{iL}, [], [], x.Quality, Path_PV_70{iL}, 1 );
end

%% Construct Figure
x.S.TraSlices       = 53;
for iL=1:2
    NativeT1w{iL}   = xASL_im_TransformData2View(IM_T1w_MNI{iL}, x);
    NativeT1w{iL}   = xASL_im_CreateVisualFig( x,  IM_T1w_MNI{iL}, [], [0.85], [], x.S.gray);
    Im50{iL}        = xASL_im_CreateVisualFig( x, {Path_T1w_MNI{iL} Path_PV_50{iL}}, [], [0.65 1], [], {x.S.gray x.S.red});
    Im70{iL}        = xASL_im_CreateVisualFig( x, {Path_T1w_MNI{iL} Path_PV_70{iL}}, [], [0.65 1], [], {x.S.gray x.S.red});
end

Column1  = [NativeT1w{1};NativeT1w{2}];
Column2  = [Im50{1};Im50{2}];
Column3  = [Im70{1};Im70{2}];
figure(1);imshow([Column1 Column2 Column3],'InitialMagnification',250,'border', 'tight')
