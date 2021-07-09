function xASL_wrp_LST_T1w_LesionFilling_WMH(x, rWMHPath)
%xASL_wrp_LST_T1w_LesionFilling_WMH Submodule of ExploreASL Structural Module, that performs lesion filling on T1w based on WMH
% segmented on the FLAIR
%
% FORMAT: xASL_wrp_LST_T1w_LesionFilling_WMH(x)
%
% INPUT:
%   x 	    - structure containing fields with all information required to run this submodule (REQUIRED)
%   x.P     - paths with NIfTIs for which this function should be applied to (REQUIRED)
%   rWMHPath - path of the WMH segmentation (either performed by LST or a copy of x.P.Path_WMH_SEGM) (REQUIRED)
%
% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This submodule runs the LST WMH-based T1w lesion filling, which should improve the registration & segmentation
% of the T1w by e.g. CAT12/SPM12. The WMH can be either segmented in the previous submodule by LST LGA/LPGA or provided externally.
% Before lesion filling, we clean up the WMH segmentation, to make the lesion filling a bit more conservative. Sometimes the WMH
% segmentation oversegments inside the GM (as there can be hyperintensities on the FLAIR) & we don't want to lesion-fill these
% on the T1w (which would turn their intensities in intensities similar to WM, leading to misclassifications by the T1w segmentation).
% Note that this is submodule only performs the lesion filling, and the clean up is also performed for the lesion filling only.
% A more thorough WMH clean up (for e.g. WMH volumetrics) is performed later in the Structural module, using also the results from the
% T1w segmentation.
%
% Note when changing the lesion filling here, LST lesion filling expects a probability map, doesnt work nicely with binary mask
% This function runs the following steps:
%
% 1. File management
% 2. Clean up the WMH segmentation used for lesion filling
% 3. Run lesion filling
% 4. Correction of too much/erronous lesion filling
% 5. File management
%
% EXAMPLE: xASL_wrp_LST_T1w_LesionFilling_WMH(x, rWMHPath);
%
% REFERENCE:
% Chard DT, Jackson JS, Miller DH, Wheeler-Kingshott CAM. Reducing the impact of white matter lesions on automated measures of brain gray and white matter volumes. J Magn Reson Imaging. 2010;32(1):223?228.
% Battaglini M, Jenkinson M, De Stefano N. Evaluating and reducing the impact of white matter lesions on brain volume measurements. Hum Brain Mapp. 2012;33(9):2062?2071.
% Pareto D, Sastre-Garriga J, Aymerich FX, et al. Lesion filling effect in regional brain volume estimations: a study in multiple sclerosis patients with low lesion load. Neuroradiology. Neuroradiology; 2016;58(5):467?474http://dx.doi.org/10.1007/s00234-016-1654-5.
% __________________________________
% Copyright 2015-2020 ExploreASL

if nargin < 2 || isempty(rWMHPath)
	error('Requires at least 2 input arguments');
elseif ~xASL_exist(rWMHPath, 'file')
    warning([rWMHPath ' missing, skipping']);
    return;
end


%% ----------------------------------------------------------------------------------
%% 1) File management
[Fpath, Ffile] = fileparts(rWMHPath);
T1_filledName = fullfile(Fpath, ['T1_filled_' Ffile(6:end) '.nii']);

fprintf('\n%s\n','----------------------------------------');
fprintf('%s\n','Removing segmented WMH from T1w, and fill lesions with values interpolated from neighborhood');


%% ----------------------------------------------------------------------------------
%% 2) Clean up the WMH segmentation used for lesion filling
xASL_im_CleanupWMHnoise(rWMHPath, rWMHPath, 200, 0.5);
% cutoff of 200 mm^3 lesion volume & pWMH>50%
% This makes sure that we only fill significant lesions


%% ----------------------------------------------------------------------------------
%% 3) Run lesion filling
if ~(xASL_stat_SumNan(xASL_stat_SumNan(xASL_stat_SumNan(xASL_io_Nifti2Im(rWMHPath))))>0)
    fprintf('No WMH lesions found, skipping T1 lesion filling\n');
    xASL_Copy(x.P.Path_T1, T1_filledName, 1);
else
    matlabbatch{1}.spm.tools.LST.filling.data = {x.P.Path_T1};
    matlabbatch{1}.spm.tools.LST.filling.data_plm = {rWMHPath};
    matlabbatch{1}.spm.tools.LST.filling.html_report = 0; % saves time
    spm_jobman('run',matlabbatch);
    close all

    %% ----------------------------------------------------------------------------------
    %% 4) Correction of too much/erronous lesion filling
    % LST lesion filling can create artifacts, which we try to remove here
    % Note that this part assumes a T1w contrast!
    T1w = xASL_io_Nifti2Im(x.P.Path_T1);
    T1wFilled = xASL_io_Nifti2Im(T1_filledName);

    % we assume that lesion filling should increase the intensity, i.e.
    % correcting the WM lesion hypointensity in the T1w to the higher WM
    % intensity. Everywhere the intensity is reduced, this is erroneous.
    if length(size(T1wFilled))~=length(size(T1w)) || ~min(size(T1wFilled)==size(T1w))
        warning('Original & filled T1 differ in size, something going wrong');
    end

    T1wFilled(T1wFilled<T1w) = T1w(T1wFilled<T1w);
    xASL_io_SaveNifti(T1_filledName, T1_filledName, T1wFilled, [], 0);
end


%% ----------------------------------------------------------------------------------
%% 5) File management
if x.settings.DELETETEMP
    xASL_delete(rWMHPath);
    xASL_adm_DeleteFileList(Fpath, '^LST_.*FLAIR\.mat$', false, [0 Inf]); % LST mat-file
end

% Rename lesion filled T1w
% If both lesion-filled & original T1 exist, rename original T1 to T1w_ORI,
% & rename the lesion-filled T1 to T1:
if  xASL_exist(x.P.Path_T1, 'file') && xASL_exist(T1_filledName, 'file')
    if  xASL_exist(x.P.Path_T1, 'file') && ~xASL_exist(x.P.Path_T1_ORI, 'file')
        % If T1w backup wasn't created through biasfield correction, create it now
        xASL_Copy(x.P.Path_T1, x.P.Path_T1_ORI);
    end
    xASL_Move(T1_filledName, x.P.Path_T1, true);
end


end
