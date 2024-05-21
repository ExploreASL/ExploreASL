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
%
%   1. Create dummy upsampled ASL scan, for registration
%   2. Reslice pGM & pWM to hi-res ASL
%   3. Estimate effective spatial resolution of ASL
%   4. Smooth pGM & pWM to this spatial resolution
%   5. Move smoothed tissue posteriors to MNI space
%
% else: run step 3 only, which will use the effective spatial resolution
% that is default for the respective sequence:
% 2D EPI: [1 1 1] * VoxelSize
% 3D GRASE: [1.1 1.1 1.38] * VoxelSize
% 3D spiral: [4.3 4.4 10.1] * VoxelSize (assuming GE uses the upsampled 2x2x4 mm 
% & run steps 1&2, but in native space these entail presmoothing &
% downsampling.
%
% EXAMPLE: xASL_wrp_PreparePV(x);
% __________________________________
% Copyright (C) 2015-2020 ExploreASL




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
if ~xASL_exist(x.P.Path_despiked_ASL4D, 'file')
    x.P.Path_despiked_ASL4D = x.P.Path_ASL4D;
end

% If the resolution of T1 partial volume maps & CBF are already the same,
% copy the PVs
if ~xASL_im_CompareNIfTIResolutionXYZ(x.P.Path_c1T1, x.P.Path_c2T1)
    warning('pGM & pWM dont have same resolution!');
else
    if xASL_im_CompareNIfTIResolutionXYZ(x.P.Path_c1T1, x.P.Path_despiked_ASL4D) && xASL_im_CompareNIfTIResolutionXYZ(x.P.Path_c2T1, x.P.Path_despiked_ASL4D)
        % copy the PV maps rather than smooth/downsample them to the ASL
        % resolution, if the resolutions are identical
        xASL_Copy(x.P.Path_c1T1, x.P.Path_PVgm, true);
        xASL_Copy(x.P.Path_c2T1, x.P.Path_PVwm, true);
        if xASL_exist(x.P.Path_c3T1)
            xASL_Copy(x.P.Path_c3T1, x.P.Path_PVcsf, true);
		end
		
		if xASL_exist(x.P.Path_WMH_SEGM, 'file')
            xASL_Copy(x.P.Path_WMH_SEGM, x.P.Path_PVwmh, true);
		end
		
        % also do this in standard space
        xASL_Copy(x.P.Pop_Path_rc1T1, x.P.Pop_Path_PV_pGM, true);
        xASL_Copy(x.P.Pop_Path_rc2T1, x.P.Pop_Path_PV_pWM, true);
		
		if xASL_exist(x.P.Path_rWMH_SEGM, 'file')
            xASL_Copy(x.P.Path_rWMH_SEGM, x.P.Pop_Path_PV_WMH_SEGM, true);
		end
        % then skip the rest of the function
        fprintf('Indentical T1 & ASL resolutions, copying PV maps\n');
        return;
    end
end
fprintf('Different T1 & ASL resolutions, adapting T1 PV maps to ASL resolution\n');

%% ------------------------------------------------------------------------------------------
%% A. STANDARD SPACE
%% ------------------------------------------------------------------------------------------

if bStandardSpace
    %% ------------------------------------------------------------------------------------------
    %% A1. Create dummy upsampled ASL scan, for registration
    % Assumes isotropic 1.5 mm resolution    
	xASL_im_Upsample(x.P.Path_PWI, x.P.Path_rPWI, [1.5 1.5 1.5],0,[50 50 50]); % adding 25 voxels at each side for each dimension, to avoid too small FoV ASL dim resulting in wrong smoothing
	
	%% ------------------------------------------------------------------------------------------
	%% A2. Reslice pGM & pWM to hi-res ASL
	% We just put it to the ASL grid, but we don't apply the non-linear transformation to the ASL,
	% because we then go directly to standard space and only apply the T1->MNI and not ASL->T1 nonlinear transformation
	xASL_spm_reslice([x.P.Path_rPWI ',1'], [x.P.Path_c1T1 ',1'], [],[], x.settings.Quality);
	xASL_spm_reslice([x.P.Path_rPWI ',1'], [x.P.Path_c2T1 ',1'], [],[], x.settings.Quality);
	if xASL_exist(x.P.Path_WMH_SEGM, 'file')
		xASL_spm_reslice([x.P.Path_rPWI ',1'], [x.P.Path_WMH_SEGM ',1'], [],[], x.settings.Quality);
	end
	
    %% ------------------------------------------------------------------------------------------
    %% A3. Define effective spatial resolution of ASL
    x = xASL_adm_DefineASLResolution(x);

    %% ------------------------------------------------------------------------------------------
    %% A4. Smooth pGM & pWM to this spatial resolution
	xASL_spm_smooth(x.P.Path_rc1T1, x.S.optimFWHM_mm, x.P.Path_rc1T1);
	xASL_spm_smooth(x.P.Path_rc2T1, x.S.optimFWHM_mm, x.P.Path_rc2T1);
	if xASL_exist(x.P.Path_WMH_SEGM, 'file')
		xASL_spm_smooth(x.P.Path_rWMH_SEGM, x.S.optimFWHM_mm, x.P.Path_rWMH_SEGM);
	end
	
	%% ------------------------------------------------------------------------------------------
	%% A5. Move smoothed tissue posteriors to MNI space
	if xASL_exist(x.P.Path_WMH_SEGM, 'file')
		InputList   = {x.P.Path_rc1T1,x.P.Path_rc2T1,x.P.Path_rWMH_SEGM};
		OutputList  = {x.P.Pop_Path_PV_pGM,x.P.Pop_Path_PV_pWM,x.P.Pop_Path_PV_WMH_SEGM};
	else
		InputList   = {x.P.Path_rc1T1,x.P.Path_rc2T1};
		OutputList  = {x.P.Pop_Path_PV_pGM,x.P.Pop_Path_PV_pWM};
	end
	
	xASL_spm_deformations(x,InputList,OutputList,4, [], [], x.P.Path_y_ASL );
	
    %% ------------------------------------------------------------------------------------------
    %% A6. Housekeeping
	if xASL_exist(x.P.Path_WMH_SEGM, 'file')
		List2Del = {x.P.Path_rc1T1 x.P.Path_rc2T1 x.P.Path_rPWI x.P.Path_rWMH_SEGM};
	else
		List2Del = {x.P.Path_rc1T1 x.P.Path_rc2T1 x.P.Path_rPWI};
	end
	
    if x.settings.DELETETEMP
        for iL=1:length(List2Del)
            xASL_delete(List2Del{iL});
        end
    end
else
    x = xASL_adm_DefineASLResolution(x); % use default effective spatial resolution
end

%% ------------------------------------------------------------------------------------------
%% B. NATIVE SPACE
%% ------------------------------------------------------------------------------------------

% We need to take into account the resolution difference between T1 and ASL.
% The correct downsampling is done by first presmoothing the image and then reslicing using normal SPM routines.


%% B1. Presmooth and transform GM and WM segmentations to the ASL space
xASL_im_PreSmooth(x.P.Path_PWI,x.P.Path_c1T1, x.P.Path_PVgm,x.S.optimFWHM_Res_mm,[],x.P.Path_mean_PWI_Clipped_sn_mat, 1);
xASL_im_PreSmooth(x.P.Path_PWI,x.P.Path_c2T1, x.P.Path_PVwm,x.S.optimFWHM_Res_mm,[],x.P.Path_mean_PWI_Clipped_sn_mat, 1);
xASL_spm_reslice(x.P.Path_PWI, x.P.Path_PVgm, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.settings.Quality, x.P.Path_PVgm);
xASL_spm_reslice(x.P.Path_PWI, x.P.Path_PVwm, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.settings.Quality, x.P.Path_PVwm);


%% B2. Presmooth and transform the CSF segmentation to the ASL space
if xASL_exist(x.P.Path_c3T1, 'file') % Check file existence as for older versions c3T1 was not always available
	xASL_im_PreSmooth(x.P.Path_PWI,x.P.Path_c3T1,x.P.Path_PVcsf,x.S.optimFWHM_Res_mm,[],x.P.Path_mean_PWI_Clipped_sn_mat, 1);
	xASL_spm_reslice(x.P.Path_PWI, x.P.Path_PVcsf, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.settings.Quality, x.P.Path_PVcsf);
end

%% B2. Presmooth and transform the WMH segmentation map to the ASL space
if xASL_exist(x.P.Path_WMH_SEGM, 'file') 
	xASL_im_PreSmooth(x.P.Path_PWI,x.P.Path_WMH_SEGM, x.P.Path_PVwmh,x.S.optimFWHM_Res_mm,[],x.P.Path_mean_PWI_Clipped_sn_mat, 1);
	xASL_spm_reslice(x.P.Path_PWI, x.P.Path_PVwmh, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.settings.Quality, x.P.Path_PVwmh);
end

%% B3. Presmooth and transform lesion & ROI masks to the ASL space.

% Note that we transform these binary masks linearly, to preserve the partial volumes.
% The maps have then have to be binarized again with a given threshold before using them again as masks.

% Create a list of lesions
LesionROIList = xASL_adm_GetFileList(fullfile(x.D.ROOT, x.SUBJECT), '(?i)^(Lesion|ROI)_(T1|FLAIR|T2)_\d*\.nii)', 'FPList', [0 Inf]);

if ~isempty(LesionROIList)
    fprintf('\n');
    warning('Lesion & ROI MASKS are transformed to the ASL native space as PROBABILISTIC MAPS!');
    fprintf('They need to be masked/binarized before using them\n');
end

% Go through all the lesions
for iROI = 1:length(LesionROIList)
    pathLesionROIT1space = xASL_adm_UnzipNifti(LesionROIList{iROI}); % Unzip and remove .gz extension if necessary
    [Fdir, Ffile, Fext] = xASL_fileparts(pathLesionROIT1space);

	% Create the path in ASL space
	pathLesionROIASLspace = fullfile(Fdir, x.SESSION, ['PV' Ffile Fext]);

	% Presmooth
	xASL_im_PreSmooth(x.P.Path_PWI, pathLesionROIT1space, pathLesionROIASLspace, x.S.optimFWHM_Res_mm, [], x.P.Path_mean_PWI_Clipped_sn_mat, 1);

	% Reslice
	xASL_spm_reslice(x.P.Path_PWI, pathLesionROIASLspace, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.settings.Quality, pathLesionROIASLspace);
end

end
