function xASL_wrp_FLAIR_BiasfieldCorrection(x)
%xASL_wrp_FLAIR_BiasfieldCorrection Submodule of ExploreASL Structural Module, that performs a biasfield correction on T1w
% & applies it on the FLAIR
%
% FORMAT: xASL_wrp_FLAIR_BiasfieldCorrection(x)
%
% INPUT:
%   x 	    - structure containing fields with all information required to run this submodule (REQUIRED)
%   x.P     - paths with NIfTIs for which this function should be applied to (REQUIRED)
%
% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This submodule performs a biasfield correction on T1w and applies it on FLAIR. This can be useful, when there are large lesions
% on the FLAIR that hamper capturing the biasfield nicely on the FLAIR itself. In such cases, the biasfield of the T1w might be easier to obtain
% and should be the same as the FLAIR, provided they are scanned in the same scan session (i.e.g same scanner, same coil).
% BE CAREFUL: this submodule assumes that the biasfields of the T1w and FLAIR are comparable, which is not the case when one of the two (or both) are
% already biasfield corrected.
%
% EXAMPLE: xASL_wrp_FLAIR_BiasfieldCorrection(x);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2019 ExploreASL
%
% 2019-05-02 HJM




%% ---------------------------------------------------------------------------------------------------
%% 1) Perform biasfield correction on T1w
xASL_spm_BiasfieldCorrection(x.P.Path_T1, x.SPMDIR, x.Quality);


%% ---------------------------------------------------------------------------------------------------
%% 2) Apply the same biasfield to the FLAIR image (NB: assumes that the FLAIR was acquired with the same scanner/coil, & hasn't been corrected
%     for its biasfield yet, also that the T1w hasn't been corrected yet

if xASL_exist(x.P.Path_FLAIR, 'file')
    % Create biasfield image
    StructIm_BF0 = xASL_io_Nifti2Im(x.P.Path_T1);
    StructIm_BF1 = xASL_io_Nifti2Im(x.P.Path_mT1);
    BiasfieldIM  = StructIm_BF1./StructIm_BF0;

    % Extrapolating biasfield
    fprintf('%s\n','Extrapolating biasfield:   ');
    BiasfieldIM = xASL_im_ndnanfilter(BiasfieldIM, 'gauss', [8 8 8]);
    xASL_io_SaveNifti(x.P.Path_mT1, x.P.Path_BiasField_T1, BiasfieldIM, [], false);

    % First resample to FLAIR space
    xASL_spm_reslice(x.P.Path_FLAIR, x.P.Path_BiasField_T1, [], [], x.Quality, x.P.Path_BiasField_FLAIR);

    % Open FLAIR & biasfield & multiply
    IM1 = xASL_io_Nifti2Im(x.P.Path_FLAIR);
    BF1 = xASL_io_Nifti2Im(x.P.Path_BiasField_FLAIR);
    IM2 = IM1.*BF1;
    xASL_io_SaveNifti(x.P.Path_FLAIR, x.P.Path_mFLAIR, IM2, 16, false);

    xASL_Move(x.P.Path_FLAIR, x.P.Path_FLAIR_ORI);
    xASL_Move(x.P.Path_mFLAIR, x.P.Path_FLAIR,1);
    xASL_delete(x.P.Path_BiasField_FLAIR);
end


%% ---------------------------------------------------------------------------------------------------
%% 3) File management
xASL_Move(x.P.Path_T1, x.P.Path_T1_ORI);
xASL_Move(x.P.Path_mT1,x.P.Path_T1,1);
xASL_delete(x.P.Path_BiasField_T1);

xASL_io_SaveNifti(x.P.Path_T1, x.P.Path_T1, xASL_io_Nifti2Im(x.P.Path_T1), 16, false); % convert to 16 bit



end
