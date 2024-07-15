function xASL_wrp_LinearReg_T1w2MNI(x, bAutoACPC)
%xASL_wrp_LinearReg_T1w2MNI Submodule of ExploreASL Structural Module, that aligns T1w with MNI
%
% FORMAT: xASL_wrp_LinearReg_T1w2MNI(x[, bAutoACPC])
%
% INPUT:
%   x 	    - structure containing fields with all information required to run this submodule (REQUIRED)
%   x.P     - paths with NIfTIs for which this function should be applied to (REQUIRED)
%   bAutoACPC - whether center of mass alignment should be performed before SPM registration (OPTIONAL, DEFAULT = true)
%
% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This submodule registers T1w linearly to the center of MNI space, a.k.a. ACPC alignment
% The same transformation is applied to all other related scans (ASL4D, M0, FLAIR, etc.)
% This facilitates MNI-based algorithms (e.g. SPM-based segmentation), and allows for visual QC with all images
% roughly in the same space. This submodule first clips high values that can bias the registration algorithm, then performs
% a center of mass-based ACPC alignment, and then several iterations of SPM coregistration.
% Assuming that this submodule is run at the start of ExploreASL, all NIfTI orientation matrices are restored before running the registration.
%
% This function performs the following steps:
% 1. Obtain path lists
% 2. Issue warnings if any of the orientation matrices were changed (e.g. in an earlier processing run)
% 3. Perform the registration
%
% EXAMPLE: xASL_wrp_LinearReg_T1w2MNI(x);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2024 ExploreASL

% Input check
if nargin<2 || isempty(bAutoACPC)
    bAutoACPC = true;
end

% Check paths
necessaryPaths = isfield(x.D,'SPMDIR') && isfield(x.P,'Path_FLAIR') && isfield(x.P,'Path_WMH_SEGM') && ...
                 isfield(x.P,'Path_ASL4D') && isfield(x.P,'Path_M0') && isfield(x.P,'Path_ASL4D_RevPE');
if ~necessaryPaths
    warning('Seemingly you are using xASL_wrp_LinearReg_T1w2MNI without defining all necessary paths...');
end


%% ---------------------------------------------------------------------------------------------------
%% 1. Obtain path lists
refPath = fullfile(x.D.SPMDIR, 'toolbox', 'OldNorm', 'T1.nii'); % = SPM8 T1 template, to make sure it is always the same
% This T1 reference image, is a blurred one. In SPM12 there are new reference images, but less blurry.
% This one empirically works fine
OtherList{1,1} = x.P.Path_FLAIR;
OtherList{end+1,1} = x.P.Path_WMH_SEGM;
OtherList{end+1,1} = x.P.Path_T1c;
OtherList{end+1,1} = x.P.Path_T2;

% Add lesion masks to the registration list
Lesion_ROI_list = xASL_adm_GetFileList(x.dir.SUBJECTDIR, ['(?i)^(Lesion|ROI)_(' x.P.STRUCT '|' x.P.FLAIR ')_\d*\.nii$'], 'FPList', [0 Inf]);
for iS=1:length(Lesion_ROI_list)
    OtherList{end+1,1} = Lesion_ROI_list{iS};
end
    
% Check for other ScanTypes that need to be in alignment of the T1w
% Here, we assume that all NIfTIs from all other scan types (e.g., ASL, fMRI, DTI, etc)
% need to remain aligned with the T1w
otherFolders = xASL_adm_GetFileList(x.dir.SUBJECTDIR, '^.*', 'FPList', [0 Inf],true); % true for folders
for iOther=1:length(otherFolders)
    DirOther = otherFolders{iOther};
    if exist(DirOther, 'dir')
        Path_Other = xASL_adm_GetFileList(DirOther, '^(?!y_).*\.nii$', 'FPList', [0 Inf]); % any DWI or func NIfTI (quick & dirty)
        for iScan=1:length(Path_Other)
            OtherList{end+1,1} = Path_Other{iScan};
        end
    end
end


%% ---------------------------------------------------------------------------------------------------
%% 2. Issue warnings if any of the orientation matrices were changed (e.g. in an earlier processing run)
% A potential error can occur if for whatever reason, in a previous rerun some where realigned (e.g., to MNI) and others were not.
% So we will now issue a warning if any of the orientation matrices were changed (e.g. in an earlier processing run).

% Note that this can be normal for FLAIR and its derivative WMH_SEGM, as these are first resampled to T1w.nii
% For T1w.nii, this can also be normal if the T1w was already realigned to MNI in a previous run.
% So we check this for ASL only, by comparing the nii.mat to nii.mat0. Either all were realigned, or none. If only some were realigned, we issue a warning.

% Get all ASL NIfTIs
checkList = [];
for iOther=1:length(OtherList)
    if xASL_exist(OtherList{iOther})
        [~, folderName] = fileparts(fileparts(OtherList{iOther}));
        if ~isempty(regexp(folderName, '^ASL_\d*$'))
            checkList{end+1} = OtherList{iOther};
        end
    end
end

% Check orientation differences for all ASL NIfTIs
% between nii.mat0 (original orientation from the scanner)
% and nii.mat (potential new orientation after registration)
% If no registration happened, these are equal
for iCheck=1:length(checkList)
    nii = xASL_io_ReadNifti(checkList{iCheck});
    matEqual(iCheck) = min(min(nii.mat==nii.mat0));
end

% Now we issue a warning if some NIfTIs had unequal orientation matrices
% — i.e., they were previously realigned —
% and others did not. In this case, matEqual has both 1s and 0s.
if length(unique(matEqual))>1
    warning('Some ASL NIfTIs were previously realigned whereas others were not!!!');
end


%% ---------------------------------------------------------------------------------------------------
%% 3. Perform the registration
%  this will be estimated on the T1w & applied to all NIfTIs
xASL_im_ClipExtremes(x.P.Path_T1, 0.999, 0, [], 1); % First we clip high vascular intensities & normalize to 4096, for more stable image contrast

if bAutoACPC % Then start with center of mass detection & realign with this
    xASL_im_CenterOfMass(x.P.Path_T1, OtherList, 15); % 15 mm minimal offset == always apply realignment if slightly off
end
% Finally, run the SPM coregistration
xASL_spm_coreg(refPath, x.P.Path_T1, OtherList, x);


end