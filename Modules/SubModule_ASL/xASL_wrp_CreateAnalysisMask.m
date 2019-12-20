function xASL_wrp_CreateAnalysisMask(x)
%xASL_wrp_CreateAnalysisMask Create analysis mask from a combination of FoV
% & removal of high and negative intravascular ASL voxels
%
% FORMAT: xASL_wrp_CreateAnalysisMask(x)
%
% INPUT:
%   x - struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT: n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function creates an analysis mask with the following steps:
%              0) Create FoV mask (native & MNI spaces)
%              1) Detect negative vascular signal (native & MNI spaces)
%              2) Detect peak vascular signal (native & MNI spaces)
%              3) Brainmasking & FoV-masking (A) native & B) MNI spaces)
%                 - Add WM vascular parts back to the mask (defined as pWM>0.8) & remove extracranial signal
%                   In the WM, negative or peak signal is more expected from
%                   noise rather than from intra-vascular signal, not many
%                   big vessels exist in the WM
%              4) Save vascular masks
%              5) Create susceptibility mask (standard space only)
%                 Here, we combine manually segmented susceptibility artifact regions in which
%                 a population-based susceptibility probability map is created
%                 This map is combined (i.e. taking the product) with the mean control & PWI
%                 intensity distribution in these regions. This product
%                 is thresholded with the average of the 75th percentile &
%                 15% of the intensity (for a bit more robustness against individual variability in sinus sizes).
%              6) Create standard space CBF_masked image to visualize masking effect
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_wrp_CreateAnalysisMask(x);
% __________________________________
% Copyright 2015-2019 ExploreASL
%
%% 0) Create native space FoV mask
% Use either original or motion estimated ASL4D
% Use despiked ASL only if spikes were detected and new file has been created
% Otherwise, despiked_raw_asl = same as original file
if ~xASL_exist(x.P.Path_despiked_ASL4D,'file')
    x.P.Path_despiked_ASL4D = x.P.Path_ASL4D;
end

FoVim = xASL_io_Nifti2Im(x.P.Path_CBF);
FoVim(:) = 1;
xASL_io_SaveNifti(x.P.Path_CBF, x.P.Path_FoV, FoVim, 8, false);
xASL_spm_deformations(x, x.P.Path_FoV, x.P.Pop_Path_FoV, 0, [], x.P.Path_mean_PWI_Clipped_sn_mat, x.P.Path_y_ASL);

%% Deal with different readouts
DoSusceptibility = true;
if strcmp(lower(x.Sequence),'2d_epi')
    Path_Template = fullfile(x.D.MapsDir,'Templates','Susceptibility_pSignal_2D_EPI.nii');
    ClipThresholdValue = 3; % 3 MAD above median
elseif strcmp(lower(x.Sequence),'3d_grase')
    Path_Template = fullfile(x.D.MapsDir,'Templates','Susceptibility_pSignal_3D_GRASE.nii');
    ClipThresholdValue = 3; % 3 MAD above median
else % for 3D spiral: don't mask susceptibility artifacts
    DoSusceptibility = false;
    ClipThresholdValue = 5; % more homogeneous image
end

%% 1) Negative vascular signal
NegativeMaskNative = xASL_im_MaskNegativeVascularSignal(x, 1); % native space
NegativeMaskMNI = xASL_im_MaskNegativeVascularSignal(x, 2); % standard space
NegativeMaskNative = xASL_im_DilateErodeFull(NegativeMaskNative,'dilate',xASL_im_DilateErodeSphere(1));
NegativeMaskMNI = xASL_im_DilateErodeFull(NegativeMaskMNI,'dilate',xASL_im_DilateErodeSphere(1));

%% 2) Detect peak vascular signal
if xASL_exist(x.P.Path_rM0,'file')
    PositiveMaskNative = xASL_im_MaskPeakVascularSignal(x.P.Path_CBF, x.P.Path_rM0, [], ClipThresholdValue);
    PositiveMaskMNI = xASL_im_MaskPeakVascularSignal(x.P.Pop_Path_qCBF, x.P.Pop_Path_M0, [], ClipThresholdValue);
else
    PositiveMaskNative = xASL_im_MaskPeakVascularSignal(x.P.Path_CBF, [], ClipThresholdValue); % no M0
    PositiveMaskMNI = xASL_im_MaskPeakVascularSignal(x.P.Pop_Path_qCBF, [], ClipThresholdValue); % no M0
end
PositiveMaskMNI = xASL_im_DilateErodeFull(PositiveMaskMNI,'dilate',xASL_im_DilateErodeSphere(1));

%% 3A) Brainmasking & FoV-masking native space
% Use previously created smoothed PGM images
pGM = xASL_io_Nifti2Im(x.P.Path_PVgm);
pWM = xASL_io_Nifti2Im(x.P.Path_PVwm);
pCSF = xASL_io_Nifti2Im(x.P.Path_PVcsf);

MaskVascularNative = ~NegativeMaskNative & ~PositiveMaskNative;
BrainMask = (pGM+pWM)>0.5;
MaskVascularNative(~BrainMask) = 0; % Remove extracranial (same setting as in ROI module)
MaskVascularNative(pWM>0.8) = 1; % Remove WM vascular spots
MaskVascularNative(pCSF>0.8) = 1; % Remove WM vascular spots

%% 3B) Brainmasking & FoV-masking standard space
pGM = xASL_io_Nifti2Im(x.P.Pop_Path_rc1T1);
pWM = xASL_io_Nifti2Im(x.P.Pop_Path_rc2T1);
pCSF = xASL_io_Nifti2Im(x.P.Pop_Path_rc3T1);
FoVim = xASL_io_Nifti2Im(x.P.Pop_Path_FoV);

MaskVascularMNI = ~NegativeMaskMNI & ~PositiveMaskMNI;
BrainMask = (pGM+pWM)>0.5 & FoVim; % -> same setting as used in ROI analysis
MaskVascularMNI(~BrainMask) = 0; % Remove extracranial & FoVim
MaskVascularMNI(pWM>0.9) = 1; % Remove WM vascular spots
MaskVascularMNI(pCSF>0.9) = 1; % Remove CSF vascular spots

%% 4) Save vascular masks
xASL_io_SaveNifti(x.P.Path_CBF, x.P.Path_MaskVascular, MaskVascularNative, 8, false);
xASL_io_SaveNifti(x.P.Pop_Path_qCBF, x.P.Pop_Path_MaskVascular, MaskVascularMNI, 8, false);

%% 5) Create susceptibility mask in standard space
if DoSusceptibility
	if xASL_exist(x.P.Pop_Path_mean_control)
		ControlIm = xASL_io_Nifti2Im(x.P.Pop_Path_mean_control); % load control image
		ControlIm = xASL_im_ndnanfilter(ControlIm, 'gauss',[2 2 2]);
	else
		ControlIm = 1;
        warning('Please check your susceptibility mask(s), we could only create it with the PWI, no control image available');
	end
     PWIIm = xASL_io_Nifti2Im(x.P.Pop_Path_PWI); % load PWI image
     PWIIm = xASL_im_ndnanfilter(PWIIm, 'gauss',[4 4 4]);
     pTemplate = xASL_io_Nifti2Im(Path_Template); % load probability map
     MaskSuscept = pTemplate<0.95*max(pTemplate(:)); % create susceptibility mask

     MixedIm = pTemplate.^0.25.*ControlIm.*PWIIm; % combine images into single probability map
     % we want to limit the influence of the template a bit, which is why
     % we use the .^0.25
     MixedIm(MixedIm<0) = 0; % clip below zero
     MixedIm(~MaskSuscept) = NaN; % remove everything outside susceptibility mask
     MixedIm = MixedIm./max(MixedIm(:)); % scale probabilities (not required)

     ThresholdPercentile = 0.75;
     % select where in the histogram we threshold
     % with proper registration, this is more dependent on smoothness of
     % image than anything else. However, this might not work with bigger
     % or smaller sinuses, which is why we also set an intensity threshold:
     ThresholdIntensity = 0.125;
     % PM: SHOULD THIS THRESHOLD TAKE THE VASCULAR MIN & MAX THRESHOLDS INTO ACCOUNT?

     SortedInt = sort(MixedIm(isfinite(MixedIm)));
     % bugfix in case of poor registration:
     if isempty(SortedInt)
         warning('Something went wrong creating susceptibility mask, poor registration?');
         return;
     end
     % don't mask for FoV, to have the same percentile everytime
     SortedInt = SortedInt(round(ThresholdPercentile*length(SortedInt)));

     FinalThreshold = (ThresholdIntensity+SortedInt)/2; % average of 2 thresholds
     FinalMask = MixedIm>FinalThreshold;

     MaskSusceptibility = BrainMask;
     MaskSusceptibility(MaskSuscept) = FinalMask(MaskSuscept);

     xASL_io_SaveNifti(x.P.Pop_Path_PWI, x.PathPop_MaskSusceptibility, MaskSusceptibility, [], false);
end

%% 6) Create standard space masked image to visualize masking effect
MaskedCBF = xASL_io_Nifti2Im(x.P.Pop_Path_qCBF);
MaskedCBF(~MaskVascularMNI) = NaN;
if DoSusceptibility
    MaskedCBF(~MaskSusceptibility) = NaN;
end
xASL_io_SaveNifti(x.P.Pop_Path_qCBF, x.P.Pop_Path_CBF_Masked, MaskedCBF, [], false);

end
