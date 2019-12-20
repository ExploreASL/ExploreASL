function x = xASL_wrp_PreparePV(x, bStandardSpace)
%xASL_wrp_PreparePV Submodule of ExploreASL ASL Module, to prepare PV maps on ASL resolution
%
% FORMAT: xASL_wrp_PreparePV(x, bStandardSpace)
%
% INPUT:
%   x               - structure containing fields with all information required to run this submodule (REQUIRED)
%   bStandardSpace  - boolean to select the creation of PV maps in standard space (OPTIONAL, DEFAULT = false)
%
% OUTPUT:
%   x               - structure containing fields with all information required to run this submodule
% OUTPUT FILES: NIfTI containing PV maps on ASL resolution
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This submodule prepares partial volume correction (PVC) by creating correct PV maps in ASL resolution, in native space,
%              as well as in standard space if requested (to perform PVC in
%              standard space):
%
% If bStandardSpace:
%   1) Create dummy upsampled ASL scan, for registration
%   2) Reslice pGM & pWM to hi-res ASL
%   3) Estimate effective spatial resolution of ASL
%   4) Smooth pGM & pWM to this spatial resolution
%   5) Move smoothed tissue posteriors to MNI space
% else: run step 3 only, which will use the effective spatial resolution
% that is default for the respective sequence:
% 2D EPI: [1 1 1] * VoxelSize
% 3D GRASE: [1.1 1.1 1.38] * VoxelSize
% 3D spiral: [4.3 4.4 10.1] * VoxelSize (assuming GE uses the upsampled 2x2x4 mm 
% & run steps 1&2, but in native space these entail presmoothing &
% downsampling
%
% EXAMPLE: xASL_wrp_PreparePV(x);
% __________________________________
% Copyright (C) 2015-2019 ExploreASL




%% ------------------------------------------------------------------------------------------
%% 0)   Admin
if nargin<2 || isempty(bStandardSpace)
	bStandardSpace = false;
end
if ~bStandardSpace
    x.ResolutionEstimation = false; % resolution estimation is currently only programmed for standard space resolution (1.5 mm)
end

% Use either original or motion estimated ASL4D
% Use despiked ASL only if spikes were detected and new file has been created
% Otherwise, despiked_raw_asl = same as original file
if ~xASL_exist(x.P.Path_despiked_ASL4D,'file')
    x.P.Path_despiked_ASL4D = x.P.Path_ASL4D;
end

if bStandardSpace
    %% ------------------------------------------------------------------------------------------
    %% 1)   Create dummy upsampled ASL scan, for registration
    % Assumes isotropic 1.5 mm resolution    
	xASL_im_Upsample(x.P.Path_PWI, x.P.Path_rPWI, [1.5 1.5 1.5],0,[50 50 50]); % adding 25 voxels at each side for each dimension, to avoid too small FoV ASL dim resulting in wrong smoothing
	
	%% ------------------------------------------------------------------------------------------
	%% 2) Reslice pGM & pWM to hi-res ASL
	xASL_spm_reslice( [x.P.Path_rPWI ',1'], [x.P.Path_c1T1 ',1'], x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.Quality);
	xASL_spm_reslice( [x.P.Path_rPWI ',1'], [x.P.Path_c2T1 ',1'], x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.Quality);

    %% ------------------------------------------------------------------------------------------
    % 3)    Estimate effective spatial resolution of ASL
    x = xASL_wrp_ResolutionEstimation(x);

    %% ------------------------------------------------------------------------------------------
    % 4)    Smooth pGM & pWM to this spatial resolution
	xASL_spm_smooth( x.P.Path_rc1T1, x.S.optimFWHM_mm, x.P.Path_rc1T1);
	xASL_spm_smooth( x.P.Path_rc2T1, x.S.optimFWHM_mm, x.P.Path_rc2T1);
	
	%% ------------------------------------------------------------------------------------------
	% 5)    Move smoothed tissue posteriors to MNI space
	InputList   = {x.P.Path_rc1T1,x.P.Path_rc2T1};
	OutputList  = {x.P.Pop_Path_PV_pGM,x.P.Pop_Path_PV_pWM};
	
	xASL_spm_deformations(x,InputList,OutputList,4, [], [], x.P.Path_y_ASL );
	% Housekeeping
	List2Del = {x.P.Path_rc1T1 x.P.Path_rc2T1 x.P.Path_rPWI};
	
	if x.DELETETEMP
		for iL=1:length(List2Del)
			xASL_delete(List2Del{iL});
		end
    end
else
    x = xASL_wrp_ResolutionEstimation(x); % use default effective spatial resolution
end

%% ------------------------------------------------------------------------------------------
% 6     Prepare pGM and pWM in the ASL native space
% by presmoothing before downsampling to ASL space
xASL_im_PreSmooth(x.P.Path_PWI,x.P.Path_c1T1,x.P.Path_sPVgm,x.S.optimFWHM_Res_mm,[],x.P.Path_mean_PWI_Clipped_sn_mat, 1);
xASL_im_PreSmooth(x.P.Path_PWI,x.P.Path_c2T1,x.P.Path_sPVwm,x.S.optimFWHM_Res_mm,[],x.P.Path_mean_PWI_Clipped_sn_mat, 1);
xASL_spm_reslice(x.P.Path_PWI, x.P.Path_sPVgm, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.Quality, x.P.Path_PVgm);
xASL_spm_reslice(x.P.Path_PWI, x.P.Path_sPVwm, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.Quality, x.P.Path_PVwm);

if xASL_exist(x.P.Path_c3T1) % for backward compatibility, where c3T1 wasnt always available
	xASL_im_PreSmooth(x.P.Path_PWI,x.P.Path_c3T1,x.P.Path_sPVcsf,x.S.optimFWHM_Res_mm,[],x.P.Path_mean_PWI_Clipped_sn_mat, 1);
	xASL_spm_reslice(x.P.Path_PWI, x.P.Path_sPVcsf, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.Quality, x.P.Path_PVcsf);
	xASL_delete(x.P.Path_sPVcsf);
end

xASL_delete(x.P.Path_sPVgm);
xASL_delete(x.P.Path_sPVwm);

end
