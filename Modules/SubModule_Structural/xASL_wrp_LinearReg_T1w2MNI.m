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
% EXAMPLE: xASL_wrp_LinearReg_T1w2MNI(x);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2019 ExploreASL
%
% 2019-05-02 HJM

% Input check
if nargin<2 || isempty(bAutoACPC)
    bAutoACPC = true;
end

% Check paths
necessaryPaths = isfield(x,'SPMDIR') && isfield(x.P,'Path_FLAIR') && isfield(x.P,'Path_WMH_SEGM') && ...
                 isfield(x.P,'Path_ASL4D') && isfield(x.P,'Path_M0') && isfield(x.P,'Path_ASL4D_RevPE');
if ~necessaryPaths
    warning('Seemingly you are using xASL_wrp_LinearReg_T1w2MNI without defining all necessary paths...');
end

%% ---------------------------------------------------------------------------------------------------
%% 1) Restore the orientation matrix of all images, in case we perform a re-run: but only when we don't have lesion maps
Lesion_ROI_list = xASL_adm_GetFileList(x.SUBJECTDIR, ['^(Lesion|ROI)_(' x.P.STRUCT '|' x.P.FLAIR ')_\d*\.nii$'], 'FPList', [0 Inf]);

%% ---------------------------------------------------------------------------------------------------
%% 2)Obtain lists of paths

refPath = fullfile(x.SPMDIR, 'toolbox','OldNorm','T1.nii'); % = SPM8 T1 template, to make sure it is always the same
% This T1 reference image, is a blurred one. In SPM12 there are new reference images, but less blurry.
% This one empirically works fine
OtherList{1,1} = x.P.Path_FLAIR;
OtherList{end+1,1} = x.P.Path_WMH_SEGM;

% Add lesion masks to the registration list
for iS=1:length(Lesion_ROI_list)
    OtherList{end+1,1} = Lesion_ROI_list{iS};
end
% Add ASL images to the registration list
for iSess = 1:x.nSessions
    OtherList{end+1,1} = x.P.Path_ASL4D;
    OtherList{end+1,1} = x.P.Path_M0;
    OtherList{end+1,1} = x.P.Path_ASL4D_RevPE;
    
    % Check for other ScanTypes
    OtherScanTypes = {'dwi' 'func'};
    for iOther=1:length(OtherScanTypes)
        DirOther = fullfile(x.SUBJECTDIR, OtherScanTypes{iOther});
        if exist(DirOther,'dir')
            Path_Other = xASL_adm_GetFileList(DirOther, '^(?!y_).*\.nii$', 'FPList', [0 Inf]); % any DWI or func NIfTI (quick & dirty)
            for iScan=1:length(Path_Other)
                OtherList{end+1,1} = Path_Other{iScan};
            end
        end
    end
end



%% ---------------------------------------------------------------------------------------------------
%% 3)Perform the registration

xASL_im_ClipExtremes(x.P.Path_T1, 0.999, 0, [], 1); % First we clip high vascular intensities & normalize to 4096, for more stable image contrast

if bAutoACPC % Then start with center of mass detection & realign with this
    xASL_im_CenterOfMass(x.P.Path_T1, OtherList, 15); % 15 mm minimal offset == always apply realignment if slightly off
end
% Finally, run the SPM coregistration
xASL_spm_coreg(refPath, x.P.Path_T1, OtherList, x);

    
    
end

