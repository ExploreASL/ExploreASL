function xASL_wrp_CreateIndividualMask(x)
%xASL_wrp_CreateIndividualMask Create analysis mask from a combination of FoV
% & removal of high and negative intravascular ASL voxels
%
% FORMAT: xASL_wrp_CreateIndividualMask(x)
%
% INPUT:
%   x - struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT: masks in native & standard space:
% native space:
%   - x.P.Path_FoV = mask of the original ASL field-of-view (FoV)
%   - x.P.Path_MaskVascular = brainmask (pGM+pWM) with vascular voxels excluded
%   - x.P.Path_MaskSusceptibility =  brainmask (pGM+pWM) with susceptibility artifact voxels excluded
%   - x.P.Path_BrainMaskProcessing = brainmask for image processing — e.g., BASIL = (pGM+pWM+pCSF).*FoV mask
%
% standard space:
%   - x.P.Pop_Path_FoV = mask of the original ASL field-of-view (FoV) transformed to standard space
%   - x.P.Pop_Path_MaskVascular = brainmask (pGM+pWM) with vascular voxels excluded
%   - x.P.Pop_Path_MaskSusceptibility = brainmask (pGM+pWM) with susceptibility artifact voxels excluded
%   - x.P.Pop_Path_BrainMaskProcessing = brainmask for image processing — e.g., BASIL = (pGM+pWM+pCSF).*FoV mask
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function creates an analysis mask with the following steps:
%              0. Create FoV mask (native & MNI spaces)
%              1. Detect negative vascular signal (native & MNI spaces, within pGM>0.5)
%              2. Detect peak vascular signal (native & MNI spaces, within pGM==80% percentile on ASL image)
%              3A. Brainmasking & FoV-masking native space
%              3B. Brainmasking & FoV-masking standard space
%              3C. Save brain mask for image processing (e.g., BASIL)
%              4. Save vascular masks
%                 - Add WM vascular parts back to the mask (defined as pWM>0.8) & remove extracranial signal
%                   In the WM, negative or peak signal is more expected from
%                   noise rather than from intra-vascular signal, not many
%                   big vessels exist in the WM
%              5. Create susceptibility mask (standard space only)
%                 Here, we combine manually segmented susceptibility artifact regions in which
%                 a population-based susceptibility probability map is created
%                 This map is combined (i.e. taking the product) with the mean control & PWI
%                 intensity distribution in these regions. This product
%                 is thresholded with the average of the 75th percentile &
%                 15% of the intensity (for a bit more robustness against individual variability in sinus sizes).
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_wrp_CreateIndividualMask(x);
% __________________________________
% Copyright 2015-2024 ExploreASL
%
%% 0. Create native space FoV mask
% Use either original or motion estimated ASL4D
% Use despiked ASL only if spikes were detected and new file has been created
% Otherwise, despiked_raw_asl = same as original file
if ~xASL_exist(x.P.Path_despiked_ASL4D,'file')
    x.P.Path_despiked_ASL4D = x.P.Path_ASL4D;
end
if ~xASL_exist(x.P.Path_PWI, 'file')
    error([x.P.Path_PWI ' missing']);
end

FoVim = xASL_io_Nifti2Im(x.P.Path_PWI);
FoVim(:) = 1;
xASL_io_SaveNifti(x.P.Path_PWI, x.P.Path_FoV, FoVim, 8, false);

if exist(x.P.Path_mean_PWI_Clipped_sn_mat, 'file') % Backwards compatability, and also needed for the Affine+DCT co-registration of ASL-T1w
    AffineTransfPath = x.P.Path_mean_PWI_Clipped_sn_mat;
else
    AffineTransfPath = [];
end

xASL_spm_deformations(x, x.P.Path_FoV, x.P.Pop_Path_FoV, 0, [], AffineTransfPath, x.P.Path_y_ASL);

%% Deal with different readouts
switch lower(x.Q.Sequence)
    case '2d_epi'
        Path_Template = fullfile(x.D.MapsDir,'Templates','Susceptibility_pSignal_2D_EPI.nii');
        ClipThresholdValue = 3; % 3 MAD above median
        DoSusceptibility = true;
    case '3d_grase'
        Path_Template = fullfile(x.D.MapsDir,'Templates','Susceptibility_pSignal_3D_GRASE.nii');
        ClipThresholdValue = 3; % 3 MAD above median
        DoSusceptibility = true;
    case '3d_spiral'
        DoSusceptibility = false;
        ClipThresholdValue = 5; % more homogeneous image        
    otherwise
        error('Unknown ASL sequence!');
end

%% 1. Negative vascular signal
NegativeMaskNative = xASL_im_MaskNegativeVascularSignal(x, 1); % native space
NegativeMaskMNI = xASL_im_MaskNegativeVascularSignal(x, 2); % standard space
NegativeMaskNative = xASL_im_DilateErodeFull(NegativeMaskNative,'dilate',xASL_im_DilateErodeSphere(1));
NegativeMaskMNI = xASL_im_DilateErodeFull(NegativeMaskMNI,'dilate',xASL_im_DilateErodeSphere(1));

%% 2. Detect peak vascular signal
if xASL_exist(x.P.Path_rM0,'file')
    PositiveMaskNative = xASL_im_MaskPeakVascularSignal(x.P.Path_PWI, x.P.Path_rM0, [], ClipThresholdValue);
    PositiveMaskMNI = xASL_im_MaskPeakVascularSignal(x.P.Pop_Path_PWI, x.P.Pop_Path_M0, [], ClipThresholdValue);
else
    PositiveMaskNative = xASL_im_MaskPeakVascularSignal(x.P.Path_PWI, [], ClipThresholdValue); % no M0
    PositiveMaskMNI = xASL_im_MaskPeakVascularSignal(x.P.Pop_Path_PWI, [], ClipThresholdValue); % no M0
end
PositiveMaskMNI = xASL_im_DilateErodeFull(PositiveMaskMNI,'dilate',xASL_im_DilateErodeSphere(1));

%% 3A. Brainmasking & FoV-masking native space
% Use previously created smoothed PGM images
pGM = xASL_io_Nifti2Im(x.P.Path_PVgm);
pWM = xASL_io_Nifti2Im(x.P.Path_PVwm);
pCSF = xASL_io_Nifti2Im(x.P.Path_PVcsf);

MaskVascularNative = ~NegativeMaskNative & ~PositiveMaskNative;
BrainMask = (pGM+pWM)>0.4;
MaskVascularNative(~BrainMask) = 0; % Remove extracranial (same setting as in ROI module)
MaskVascularNative(pWM>0.8) = 1; % Remove WM vascular spots
MaskVascularNative(pCSF>0.8) = 1; % Remove WM vascular spots

% Obtain brain mask for image processing (e.g., BASIL)
BrainMaskProcessingNativeSpace = (pGM+pWM+pCSF)>0.08;

%% 3B. Brainmasking & FoV-masking standard space
pGM = xASL_io_Nifti2Im(x.P.Pop_Path_PV_pGM);
pWM = xASL_io_Nifti2Im(x.P.Pop_Path_PV_pWM);
pCSF = xASL_io_Nifti2Im(x.P.Pop_Path_PV_pCSF);
FoVim = xASL_io_Nifti2Im(x.P.Pop_Path_FoV);

MaskVascularMNI = ~NegativeMaskMNI & ~PositiveMaskMNI;
BrainMask = (pGM+pWM)>0.4 & FoVim; % -> same setting as used in ROI analysis
MaskVascularMNI(~BrainMask) = 0; % Remove extracranial & FoVim
MaskVascularMNI(pWM>0.9) = 1; % Remove WM vascular spots
MaskVascularMNI(pCSF>0.9) = 1; % Remove CSF vascular spots

% Obtain brain mask for image processing (e.g., BASIL)
BrainMaskProcessingStandardSpace = (pGM+pWM+pCSF)>0.08 & FoVim;

%% 3C. Save brain mask for image processing (e.g., BASIL)
% This mask can be used for fitting data, e.g., BASIL, fitting ATT, Tex,
% etc. that is not useful in extracranial voxels, and beneficial to speed
% up by only processing intracranial voxels. On the other hand, it can be
% helpful to also process some voxels outside the brain for inspecting
% artifacts, and registration errors shouldn't lead to accidentally
% excluding intracranial voxels. Therefore, we use a rather inclusive mask.

xASL_io_SaveNifti(x.P.Path_PWI, x.P.Path_BrainMaskProcessing, BrainMaskProcessingNativeSpace, 8, false);
xASL_io_SaveNifti(x.P.Pop_Path_PWI, x.P.Pop_Path_BrainMaskProcessing, BrainMaskProcessingStandardSpace, 8, false);

%% 4. Save vascular masks
xASL_io_SaveNifti(x.P.Path_PWI, x.P.Path_MaskVascular, MaskVascularNative, 8, false);
xASL_io_SaveNifti(x.P.Pop_Path_PWI, x.P.Pop_Path_MaskVascular, MaskVascularMNI, 8, false);

%% 5. Create susceptibility mask in standard space
if DoSusceptibility
	
    if xASL_exist(x.P.Pop_Path_noSmooth_M0)
        ControlIm = xASL_io_Nifti2Im(x.P.Pop_Path_noSmooth_M0); % load M0 image
    elseif xASL_exist(x.P.Pop_Path_mean_control)
		ControlIm = xASL_io_Nifti2Im(x.P.Pop_Path_mean_control); % load control image
	else
		ControlIm = 1;
        warning('Please check your susceptibility mask(s), we could only create it with the PWI, no M0 or control image available');
    end
    
     if numel(ControlIm)~=1
         ControlIm = xASL_im_ndnanfilter(ControlIm, 'gauss',[2 2 2]);
     end
    
     PWIIm = xASL_io_Nifti2Im(x.P.Pop_Path_PWI); % load PWI image
     PWIIm = xASL_im_ndnanfilter(PWIIm, 'gauss',[4 4 4]);
     pTemplate = xASL_io_Nifti2Im(Path_Template); % load probability map
     
     % Change pTemplate based on sequence (thanks to Khazar for checking
     % this for 3D GRASE in the BioFinder study)
     switch lower(x.Q.Sequence)
         case '2d_epi'
             SusceptibilityThreshold = 0.95;
         case '3d_grase'
             SusceptibilityThreshold = 0.6;
         otherwise
             warning('Unknown sequence for susceptibility thresholding, skipping');
             SusceptibilityThreshold = 1;
     end
	 MaskSuscept = pTemplate<SusceptibilityThreshold*max(pTemplate(:)); % create susceptibility mask

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

     % %% Brief explanation what these masks here mean:
     % BrainMask   => brain parenchyma (GM+WM) == 1, CSF & outside the brain == 0
     % MaskSuscept => potential susceptibility regions, where artifacts may occur == 1, otherwise == 0
     % FinalMask   => within MaskSuscept, voxels we want to include in the mask == 1, otherwise == 0

     MaskSusceptibility = BrainMask; % we start with a wholebrain parenchyme mask
     MaskSusceptibility(MaskSuscept) = FinalMask(MaskSuscept); % within potential susceptibility regions, we mask out susceptibility artifacts
     MaskSusceptibility(~BrainMask) = 0; % we ensure that voxels without susceptibility artifacts, outside the brain, are still masked out
else
    MaskSusceptibility = BrainMask;
    % for e.g. 3D spiral the susceptibility mask is equal to the brain
    % mask, effectively not masking out susceptibility artifacts
end

xASL_io_SaveNifti(x.P.Pop_Path_PWI, x.P.Pop_Path_MaskSusceptibility, MaskSusceptibility, [], false);


end