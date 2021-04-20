function [QC] = XQC_Noise_IT(NiftiPath, Pathc1T1, Pathc2T1, SubjDir, AtlasDir, x, y_file)

%% Load Settings and Images

% Load all the files
NiftiHDR = xASL_io_ReadNifti(NiftiPath);
Im = xASL_io_Nifti2Im(NiftiPath);

%Reslice GM and WM to Nifti
xASL_spm_reslice(NiftiPath,Pathc1T1 , [], [], [], fullfile(SubjDir, 'tmp_GM_mask.nii'), 0)
xASL_spm_reslice(NiftiPath,Pathc2T1 , [], [], [], fullfile(SubjDir, 'tmp_WM_mask.nii'), 0)

%Load GM and WM
pGM = xASL_io_Nifti2Im(fullfile(SubjDir, 'tmp_GM_mask.nii'));
pWM = xASL_io_Nifti2Im(fullfile(SubjDir, 'tmp_WM_mask.nii'));

%Delete temporary resliced files
xASL_delete(fullfile(SubjDir, 'tmp_GM_mask.nii'))
xASL_delete(fullfile(SubjDir, 'tmp_GM_mask.nii'))


%Manage dimensionality, average over the 4Th dimension to handle 4D data
Im = xASL_stat_MeanNan(Im, 4);

% Load all the general masks
Path_DeepWMmask = fullfile(AtlasDir,'CentralWM_QC.nii');
Path_NativeDeepWM = fullfile(SubjDir,'CentralWM_QC.nii');
Path_LRMask = fullfile(AtlasDir,'LeftRight.nii');
Path_NativeLRMask = fullfile(SubjDir,'LeftRight.nii');
%Path_MNI_DeepWM = fullfile(x.D.PopDir,['CentralWM_QC_' x.P.SubjectID '.nii']); % May be useless


%% Create basic masks
% Define a simple native space large brain mask
BrainMask = (pGM+pWM)>0.35; % bit lower threshold than normal (0.5) because of potential atrophy
pCSF = max(0, BrainMask - pGM - pWM);

% Threshold the anatomical volume (in case of resampling causing negative values)
Im(Im<eps) = eps;

% make gray/white mutually exclusive
GMmask = logical(pGM>pWM & pGM>pCSF) & BrainMask;
WMmask = logical(pWM>pGM & pWM>pCSF) & BrainMask;








%% Erode WM as proxy-background noise region
% These files need to be inversely transformed from MNI to native space.
% We reuse these masks if they already exist, which speeds up QC & keeps
% the mask identical for the T1w & FLAIR

if xASL_exist(Path_NativeDeepWM,'file') && xASL_exist(Path_NativeLRMask,'file')
    xASL_spm_reslice(NiftiPath, Path_NativeDeepWM, [], [], [], Path_NativeDeepWM, 0);
    xASL_spm_reslice(NiftiPath, Path_NativeLRMask, [], [], [], Path_NativeLRMask, 0);
    WMeroded = xASL_io_Nifti2Im(Path_NativeDeepWM);
else
    % compute inverse transformation & put the MNI ROIs in native space
    xASL_spm_deformations(x ,Path_DeepWMmask, Path_NativeDeepWM, [] , NiftiPath, [], y_file, []);
    xASL_spm_deformations(x, Path_LRMask, Path_NativeLRMask, [], NiftiPath, [], y_file, [])
    IM_DeepWM = xASL_io_Nifti2Im(Path_NativeDeepWM);

    % 2) Multiply native space DeepWM template mask by subject-wise pWM=100%, & need real data (isfinite)
    WMeroded = (pWM==1) .* IM_DeepWM .* BrainMask;
    xASL_io_SaveNifti(Path_NativeDeepWM, Path_NativeDeepWM, WMeroded, [], false);
    % save mask for later use, & for visualization in standard space
end
 


%WM region stats
WMRefData = Im;
WMRefData(~WMeroded) = NaN;
MeanWMRef = xASL_stat_MeanNan(WMRefData(:));
VoxelSize = NiftiHDR.hdr.pixdim(2:4);


% High-pass filter, ignoring NaNs (voxels outside mask)
SmoothWMRefData = xASL_im_ndnanfilter(WMRefData,'gauss',[8 8 8]./double(VoxelSize),1);
WMRefData = MeanWMRef.*WMRefData ./ SmoothWMRefData;
WMRefmask = WMRefData>0;

%% Compute ROIs Values


BrainRoi = Im(BrainMask); % Whole Brain 
GMroi = Im(GMmask); % Grey Matter
WMroi = Im(WMmask); % White Matter

MeanGM = mean(GMroi);
MeanWM = mean(WMroi);

StdWMRef = xASL_stat_StdNan(WMRefData(:));
VarWMRef = xASL_stat_VarNan(WMRefData(:));




%% Extract QC Parameters
% Mean GM / SD WMref  (higher = better)
QC.Noise.SNR_GM_Ratio = MeanGM/StdWMRef;

% abs(mean GM - mean WM)/variance WMref (higher = better)
QC.Noise.CNR_GM_WM_Ratio = abs(MeanGM-MeanWM)/VarWMRef;

% Variance within brain / variance WMref (higher = better)
QC.IT.FBER_WMref_Ratio = var(BrainRoi)/VarWMRef;

% Shannon entropy of voxel intensities proportional to maximum possible entropy for similarly sized image
% indicating ghosting and head motion-induced blurring (lower = better)
Bmax = sqrt(xASL_stat_SumNan(BrainRoi(:).^2));
QC.IT.EFC_bits = -xASL_stat_SumNan((BrainRoi(:)./Bmax) .* log((BrainRoi(:)./Bmax)));

end 


