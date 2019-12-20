%% Create sequence-based susceptibility artifact probability map

% WBmask
pGM = xASL_io_Nifti2Im('C:\BackupWork\ASL\EPAD_xASL_Paper\analysis_EPAD\Population\Templates\PV_pGM_bs-mean.nii');
pWM = xASL_io_Nifti2Im('C:\BackupWork\ASL\EPAD_xASL_Paper\analysis_EPAD\Population\Templates\PV_pWM_bs-mean.nii');

WBmask = (pGM+pWM)>0.25;

% MNI structural
MNIpath = 'C:\ExploreASL\Maps\Atlases\MNI_structural.nii';
SubCort = xASL_io_Nifti2Im(MNIpath);

% Subcortical GM structures
% 9 thalamus
% 7 putamen
% 1 caudate

SubCort = SubCort==1 | SubCort==7 | SubCort==9;
SubCort = SubCort.*pGM;

%% Define susceptibility regions
% Open EPAD T1 for visualization
clear IM
IM = xASL_io_Nifti2Im('C:\BackupWork\ASL\EPAD_xASL_Paper\analysis_EPAD\Population\Templates\T1_bs-mean_Unmasked.nii');
Mask = zeros(size(IM));
% Define mastoid air cavity
Mask([15:45 121-45:121-15],45:75,12:40) = 1;

% Frontal sinus
Mask(40:121-40,125:140,30:50) = 1;

% Ethmoid
Mask(40:121-40,70:125,25:38) = 1;

% Mask out WM
Mask(pWM>0.5) = 0;
Mask = Mask.*(1-pWM);

% Mask out GM in subcortical structures
Mask(SubCort>0.5) = 0;
Mask = Mask.*(1-SubCort);


%% Create probability map
% Open control SNR image 2D EPI EPAD
clear IM
IM{1} = xASL_io_Nifti2Im('C:\BackupWork\ASL\EPAD_xASL_Paper\analysis_EPAD\Population\Templates\SNR_control_bs-mean_Unmasked.nii');
IM{2} = xASL_io_Nifti2Im('C:\BackupWork\ASL\Harmy\analysis_Harmy\Population\Templates\SNR_control_bs-mean.nii');

for iIM=1:length(IM)
    IMwb{iIM} = IM{iIM}(WBmask);
    ThresholdIM(iIM) = median(IMwb{iIM}(:));
    pIM{iIM} = IM{iIM};
    pIM{iIM}(pIM{iIM}>ThresholdIM(iIM)) = ThresholdIM(iIM);
    pIM{iIM} = pIM{iIM}./ThresholdIM(iIM);
    pIM{iIM}(~Mask) = 1;
    pIM{iIM} = xASL_im_ndnanfilter(pIM{iIM},'gauss',[6 6 6]);  % smooth
end

%% Save probability maps
Path_EPImask = 'C:\ExploreASL\Maps\Templates\Susceptibility_pSignal_2D_EPI.nii';
Path_GRASEmask = 'C:\ExploreASL\Maps\Templates\Susceptibility_pSignal_3D_GRASE.nii';

xASL_io_SaveNifti(MNIpath, Path_EPImask, pIM{1}, [], false);
xASL_io_SaveNifti(MNIpath, Path_GRASEmask, pIM{2}, [], false);
