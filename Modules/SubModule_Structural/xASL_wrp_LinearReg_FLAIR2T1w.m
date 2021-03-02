function xASL_wrp_LinearReg_FLAIR2T1w(x, bAutoACPC)
%xASL_wrp_LinearReg_FLAIR2T1w Submodule of ExploreASL Structural Module, that aligns FLAIR with T1w
%
% FORMAT: xASL_wrp_LinearReg_FLAIR2T1w(x[, bAutoACPC])
%
% INPUT:
%   x 	    - structure containing fields with all information required to run this submodule (REQUIRED)
%   x.P     - paths with NIfTIs for which this function should be applied to (REQUIRED)
%   bAutoACPC - whether center of mass alignment should be performed before SPM registration (OPTIONAL, DEFAULT = true)
%
% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This submodule registers FLAIR linearly to the T1w
% The same transformation is applied to all other related scans (FLAIR-segmented lesions, WMH specifically or other lesions)
% This is required to enable the application of T1w derivatives (e.g. transformations to standard space, tissue segmentation) for FLAIR 
% and vice versa (e.g. WMH lesion-filling).
%
% EXAMPLE: xASL_wrp_LinearReg_FLAIR2T1w(x);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2019 ExploreASL
%
% 2019-05-02 HJM


if nargin<2 || isempty(bAutoACPC)
    bAutoACPC = true;
end

if ~xASL_exist(x.P.Path_FLAIR, 'file')
	return;
end

fprintf('\n%s\n','---------------------------------');
fprintf('%s\n','FLAIR.nii detected, processing...');


%% ---------------------------------------------------------------------------------------------------
%% 1)Obtain lists of paths

OtherList = {x.P.Path_WMH_SEGM};

Lesion_FLAIR_list = xASL_adm_GetFileList(x.SUBJECTDIR, ['^Lesion_' x.P.FLAIR '_\d*\.(nii|nii\.gz)$'], 'FPList', [0 Inf]);
for iS=1:length(Lesion_FLAIR_list)
    OtherList{end+1,1} = Lesion_FLAIR_list{iS};
end



%% ---------------------------------------------------------------------------------------------------
%% 2)Perform the registration

xASL_im_ClipExtremes(x.P.Path_FLAIR, 0.999999, 0, [], 1); % First we clip high vascular intensities & normalize to 4096, for more stable image contrast

if bAutoACPC % Then start with center of mass detection & realign with this
    xASL_im_CenterOfMass(x.P.Path_FLAIR, OtherList);
end
% Finally, run the SPM coregistration
xASL_spm_coreg(x.P.Path_T1, x.P.Path_FLAIR, OtherList, x);

    
end

