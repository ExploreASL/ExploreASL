function xASL_wrp_LinearReg_Others2T1w(x, bAutoACPC)
%xASL_wrp_LinearReg_Others2T1w Submodule of ExploreASL Structural Module, that aligns T2 and T1c with T1w
%
% FORMAT: xASL_wrp_LinearReg_Others2T1w(x[, bAutoACPC])
%
% INPUT:
%   x 	    - structure containing fields with all information required to run this submodule (REQUIRED)
%   x.P     - paths with NIfTIs for which this function should be applied to (REQUIRED)
%   bAutoACPC - whether center of mass alignment should be performed before SPM registration (OPTIONAL, DEFAULT = true)
%
% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This submodule registers T1c and T2 linearly to the T1w
%
% EXAMPLE: xASL_wrp_LinearReg_Others2T1w(x);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2021 ExploreASL

if nargin<2 || isempty(bAutoACPC)
    bAutoACPC = true;
end

if (~xASL_exist(x.P.Path_T1c, 'file')) && (~xASL_exist(x.P.Path_T2, 'file'))
	return;
end

fprintf('\n%s\n','----------------------------------------');
fprintf('%s\n','T1c.nii or T2.nii detected, processing...');


%% ---------------------------------------------------------------------------------------------------
%% 1)Perform the registration

listPaths = {x.P.Path_T2, x.P.Path_T1c};
for iFile = 1:length(listPaths)
	if xASL_exist(listPaths{iFile},'file')
		xASL_im_ClipExtremes(listPaths{iFile}, 0.999999, 0, [], 1); % First we clip high vascular intensities & normalize to 4096, for more stable image contrast
		
		if bAutoACPC % Then start with center of mass detection & realign with this
			xASL_im_CenterOfMass(listPaths{iFile}, {});
		end
		% Finally, run the SPM coregistration
		xASL_spm_coreg(x.P.Path_T1, listPaths{iFile}, {}, x);
	end
end
    
end
