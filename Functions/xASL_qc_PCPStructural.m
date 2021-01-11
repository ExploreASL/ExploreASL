function [anatQA] = xASL_qc_PCPStructural(PathT1, Pathc1T1, Pathc2T1, x, PopPathT1)
% Computes anatomical QC parameters
%
% FORMAT: [anatQA] = xASL_qc_PCPStructural(PathT1, Pathc1T1, Pathc2T1, x, PopPathT1)
%
% INPUT:
%   PathT1     - path to the anatomical image (REQUIRED)
%   Pathc1T1   - path to the GM segmentation image (REQUIRED)
%   Pathc2T1   - path to the WM segmentation image (REQUIRED)
%   x          - structure containing fields with all information required to run this submodule (REQUIRED)
%   PopPathT1  - path to the anatomical image in standard space (REQUIRED)
%
% OUTPUT:
%   anatQA     - struct containing the QC parameters
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function computes several anatomical QC parameters as proposed in SPM Univariate Plus:
%
%              - WM_ref_vol_mL    - volume of the WM reference region (mL)
%              - WMref_vol_Perc   - same but as percentage of total WM volume
%              - SNR_GM           - GM signal-to-Noise Ratio (SNR), ie the mean intensity within GM divided
%                                   by SD of WM reference region. Higher = better.
%              - CNR_GM_WM        - GM-WM Contrast-to-Noise Ratio (CNR), i.e. the mean of GM - mean of WM
%                                   divided by the SD of the WM reference region. Higher = better.
%              - FBER_WMref_Ratio - Foreground to Background Energy Ratio (FBER), i.e. the variance of voxels within the brain (in pGM+pWM mask)
%                                   divided by the variance of voxels in the WM reference region. Higher = better.
%              - EFC_bits         - Shannon Entropy Focus Criterion (EFC), i.e. the entropy of voxel intensities proportional to the maximum
%                                   possibly entropy for a similarly sized image. Indicates ghosting and head motion-induced blurring. Lower = better.
%              - Mean_AI_Perc     - mean relative voxel-wise absolute Asymmetry Index (AI) within the brain (pGM+pWM mask) (%)
%              - SD               - same but SD (%)
%
% REFERENCES:
%              Preprocessed Connectome Project Quality Assurance Protocol (QAP):
%              http://preprocessed-connectomes-project.org/quality-assessment-protocol/
%              http://ieeexplore.ieee.org/document/650886/
%
% EXAMPLE: anatQA = xASL_qc_PCPStructural('/analysis/Subj001/T1.nii', '/analysis/Subj001/c1T1.nii', '/analysis/Subj001/c2T1.nii', x, 'analysis/Population/T1.nii');
% __________________________________
% Copyright (C) 2015-2019 ExploreASL


%% Admin
if nargin<5
	error('Functions requires all five input parameters.');
end

if isempty(PathT1) || isempty(Pathc1T1) || isempty(Pathc2T1) || isempty(PopPathT1)
	error('All paths need to be defined.');
end

% Load all the files
T1nii = xASL_io_ReadNifti(PathT1);
T1im = xASL_io_Nifti2Im(PathT1);
pGM = xASL_io_Nifti2Im(Pathc1T1);
pWM = xASL_io_Nifti2Im(Pathc2T1);

% Load all the general masks
Path_DeepWMmask = fullfile(x.D.MapsSPMmodifiedDir,'CentralWM_QC.nii');
Path_NativeDeepWM = fullfile(x.SUBJECTDIR,'CentralWM_QC.nii');
Path_LRMask = fullfile(x.D.MapsSPMmodifiedDir,'LeftRight.nii');
Path_NativeLRMask = fullfile(x.SUBJECTDIR,'LeftRight.nii');
Path_MNI_DeepWM = fullfile(x.D.PopDir,['CentralWM_QC_' x.P.SubjectID '.nii']);

fprintf('%s\n','Computing anatomical QC parameters from the Preprocessed Connectome Project (PCP) QA Protocol');

%% Create basic masks
% Define a simple native space large brain mask
BrainMask = (pGM+pWM)>0.35; % bit lower threshold than normal (0.5) because of potential atrophy
pCSF = max(0, BrainMask - pGM - pWM);

% Threshold the anatomical volume (in case of resampling causing negative values)
T1im(T1im<eps) = eps;

% make gray/white mutually exclusive
GMmask = logical(pGM>pWM & pGM>pCSF) & BrainMask;
WMmask = logical(pWM>pGM & pWM>pCSF) & BrainMask;

%% Erode WM as proxy-background noise region
% These files need to be inversely transformed from MNI to native space.
% We reuse these masks if they already exist, which speeds up QC & keeps
% the mask identical for the T1w & FLAIR

if xASL_exist(Path_NativeDeepWM,'file') && xASL_exist(Path_NativeLRMask,'file')
    xASL_spm_reslice(PathT1, Path_NativeDeepWM, [], [], x.Quality, Path_NativeDeepWM, 0);
    xASL_spm_reslice(PathT1, Path_NativeLRMask, [], [], x.Quality, Path_NativeLRMask, 0);
    WMeroded = xASL_io_Nifti2Im(Path_NativeDeepWM);
else
    % compute inverse transformation & put the MNI ROIs in native space
    xASL_spm_deformations(x,{Path_DeepWMmask Path_LRMask},{Path_NativeDeepWM Path_NativeLRMask},0, PathT1);

    IM_DeepWM = xASL_io_Nifti2Im(Path_NativeDeepWM);

    % 2) Multiply native space DeepWM template mask by subject-wise pWM=100%, & need real data (isfinite)
    WMeroded = (pWM==1) .* IM_DeepWM .* BrainMask;
    xASL_io_SaveNifti(Path_NativeDeepWM, Path_NativeDeepWM, WMeroded, [], false);
    % save mask for later use, & for visualization in standard space
end

% Put native space deep WM ROI in standard space for visualization
% if it wasnt already there.
if ~xASL_exist(Path_MNI_DeepWM,'file')
    xASL_spm_deformations(x, Path_NativeDeepWM, Path_MNI_DeepWM, 0);
end

% Calculate WMref region data & stats
WMRefData = T1im;
WMRefData(~WMeroded) = NaN;
MeanWMRef = xASL_stat_MeanNan(WMRefData(:));
VoxelSize = T1nii.hdr.pixdim(2:4);

% 3) High-pass filter, ignoring NaNs (voxels outside mask)
SmoothWMRefData = xASL_im_ndnanfilter(WMRefData,'gauss',[8 8 8]./double(VoxelSize),1);
WMRefData = MeanWMRef.*WMRefData ./ SmoothWMRefData;
WMRefmask = WMRefData>0;
anatQA.WMref_vol_mL = sum(WMRefmask(:)).*prod(VoxelSize)/1000; % in mL

% WMvolume = sum(pWM(:)).*prod(VoxelSize)/1000; % in mL
% anatQA.WMref_vol_Perc = 100*(anatQA.WMref_vol_mL/WMvolume);

%% Load native space stuff & acquire ROI data
LeftRightIm = xASL_io_Nifti2Im(Path_NativeLRMask);

T1roi = T1im(isfinite(T1im)); % complete T1w
GMroi = T1im(GMmask & isfinite(T1im)); % within GM mask only
WMroi = T1im(WMmask & isfinite(T1im)); % within WM mask only

LeftMask = WMRefmask & LeftRightIm==1; % Left WMRef (i.e. deep WM) mask
RightMask = WMRefmask & LeftRightIm==2; % Left WMRef (i.e. deep WM) mask
WMRefData_Left = WMRefData(LeftMask);
WMRefData_Right = WMRefData(RightMask);
WeightLeft = sum(LeftMask(:))  / (sum(LeftMask(:)) + sum(RightMask(:)));
WeightRight = sum(RightMask(:)) / (sum(LeftMask(:)) + sum(RightMask(:)));

MeanGM = mean(GMroi);
MeanWM = mean(WMroi);
StdWMRef = std(WMRefData_Left(:))*WeightLeft + std(WMRefData_Right(:))*WeightRight;
VarWMRef = var(WMRefData_Left(:))*WeightLeft + var(WMRefData_Right(:))*WeightRight;

%% Calculate the QA parameters

anatQA.SD_WMref = StdWMRef;

% Mean GM / SD WMref  (higher = better)
anatQA.SNR_GM_Ratio = MeanGM/StdWMRef;

% abs(mean GM - mean WM)/variance WMref (higher = better)
anatQA.CNR_GM_WM_Ratio = abs(MeanGM-MeanWM)/VarWMRef;

% Variance within brain / variance WMref (higher = better)
anatQA.FBER_WMref_Ratio = var(T1roi)/VarWMRef;

% Shannon entropy of voxel intensities proportional to maximum possible entropy for similarly sized image
% indicating ghosting and head motion-induced blurring (lower = better)
Bmax = sqrt(xASL_stat_SumNan(T1roi(:).^2));
anatQA.EFC_bits = -xASL_stat_SumNan((T1roi(:)./Bmax) .* log((T1roi(:)./Bmax)));


%% Compute asymmetry in standard space (avoids registration issues)
T1im = xASL_io_Nifti2Im(PopPathT1);
pGM = xASL_io_Nifti2Im(x.P.Pop_Path_rc1T1);
pWM = xASL_io_Nifti2Im(x.P.Pop_Path_rc2T1);

BrainMask = (pGM+pWM)>0.5;
T1im(~BrainMask) = NaN;
T1im = xASL_im_ndnanfilter(T1im,'gauss',[4 4 4]./[1.5 1.5 1.5],1);

L = T1im;
R = xASL_im_Flip(T1im,1);
DiffIm = abs(L-R);
MeanIm = (L+R)./2;
asym = DiffIm./MeanIm;
anatQA.Mean_AI_Perc = xASL_stat_MeanNan(asym(:)).*100; % in percentages
anatQA.SD_AI_Perc = xASL_stat_StdNan(asym(:)).*100; % in percentages

return
