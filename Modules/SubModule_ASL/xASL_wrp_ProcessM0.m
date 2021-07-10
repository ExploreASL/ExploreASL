function xASL_wrp_ProcessM0(x)
%xASL_wrp_ProcessM0 Submodule of ExploreASL ASL Module, for M0 image processing
%
% FORMAT: xASL_wrp_ProcessM0(x)
%
% INPUT:
%   x       - structure containing fields with all information required to run this submodule (REQUIRED)
%   x.bRegisterM02ASL - boolean specifying whether M0 is registered to
%                     mean_control image (or T1w if no control image exists)
%                     It can be useful to disable M0 registration if the
%                     ASL registration is done based on the M0, and little
%                     motion is expected between the M0 and ASL
%                     acquisition.
%                     If no separate M0 image is available, this parameter
%                     will have no effect. This option is disabled
%                     automatically for 3D spiral
%                     (OPTIONAL, DEFAULT = 0)
%                     - 0 = M0 registration disabled
%                     - 1 = M0 registration enabled (DEFAULT)
%
% OUTPUT: n/a
% OUTPUT FILES: NIfTI containing image processed M0 map in native & standard space, with and without smoothing
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This submodule performs the image processing and
%           quantification of M0 maps (if they exist), with the following steps:
%
%           1. Register M0 to mean control if it exists
%              Before registration, contrast is equalized between the
%              images & biasfields are removed
%           2. Quantify M0 (correction scale slope & incomplete T1 recovery)
%           3. Masking & smoothing of M0 image, either using:
%              A) traditional technique (very sharp masking & little smoothing)
%              B) new ExploreASL-specific technique:
%                 * extrapolating outside mask (avoiding artifacts from too
%                   much or too little masking)
%                 * smooth very extensively, to create a biasfield
%                   (increases robustness & comparison of M0 between
%                   sequences/patients)
%
%     Any M0 will be processed here. Even if part of the subjects does not
%     have an M0, since this can be later imputed, or an average population
%     M0 image could be used. Also, without background suppression and
%     without an M0, the MeanControl image is before saved as M0, and will
%     be processed here as well.
%
%     Note that any voxel-size differences between M0 and ASL are allowed
%     here: step 0B below rescales the PD inside an M0 voxel to the same as
%     the ASL resolution (assuming a voxel with half volume contains half
%     the amount of protons). The M0 is further processed in standard
%     space, and reduced to a biasfield. For the quantification in standard
%     space, the PWI and M0 are now by definition in the same space.
%     Also, the standard space M0 biasfield is resampled to the native PWI
%     space (at the end of step 3B), ensuring that both are also in the
%     same native space.
%
% EXAMPLE: xASL_wrp_ProcessM0(x);
% __________________________________
% Copyright (C) 2015-2020 ExploreASL



%% -----------------------------------------------------------------------------------------------
%% 0)   Administration
if ~xASL_exist(x.P.Path_M0,'file')
    % skip this part. This should be only the case when using a single value for M0,
    % if the mean_control is used as M0, it should be copied to be used as M0 previously
    return;
end

% Use either original or motion estimated ASL4D
% Use despiked ASL only if spikes were detected and new file has been created
% Otherwise, despiked_raw_asl = same as original file
if ~xASL_exist(x.P.Path_despiked_ASL4D,'file')
    x.P.Path_despiked_ASL4D = x.P.Path_ASL4D;
end
tempnii = xASL_io_ReadNifti(x.P.Path_despiked_ASL4D);
nVolumes = double(tempnii.hdr.dim(5));

if strcmpi(x.M0,'no_background_suppression')
    x.M0 = 'UseControlAsM0'; % backward compatibility
end
if ~isfield(x.settings,'M0_conventionalProcessing')
       x.settings.M0_conventionalProcessing   = 0;
       % by default, conventional processing is off, since our new method outperforms in most cases
elseif x.settings.M0_conventionalProcessing == 1 && strcmpi(x.Q.readoutDim,'3D')
       x.settings.M0_conventionalProcessing = 0;
       warning('M0 conventional processing disabled, since this masking does not work with 3D sequences');
end

%% 0B) Scale PD between ASL & M0 if voxel sizes differ
M0nii       = xASL_io_ReadNifti(x.P.Path_M0);
M0size      = prod(M0nii.hdr.pixdim(2:4));
ASLnii      = xASL_io_ReadNifti(x.P.Path_ASL4D);
ASLsize     = prod(ASLnii.hdr.pixdim(2:4));
M0ScaleF    = ASLsize/M0size;

% copy existing M0 for processing, in single precision
% averaging if multiple frames, as we will blur later anyway,
% we can skip the motion correction here
xASL_io_SaveNifti(x.P.Path_M0, x.P.Path_rM0, mean(xASL_io_Nifti2Im(x.P.Path_M0).*M0ScaleF,4), 32, 0);

% Note that here we created rM0, which is averaged across 4th dimension, and adapted along this function

xASL_im_CreateASLDeformationField(x); % make sure we have the deformation field in ASL resolution


%% -----------------------------------------------------------------------------------------------
%% 1)   Register M0 to mean control if it exists
% Registering x.P.Path_M0 to ASL, changing x.P.Path_M0 header only
% If there is only a single ASL PWI (e.g. GE 3D FSE), this is not performed because of
% inequality of image contrast for ASL & M0, and because they usually
% are already in decent registration.

if isfield(x, 'bRegisterM02ASL') && ~x.bRegisterM02ASL
    fprintf('M0 registration (to ASL or T1w) is skipped upon request\n');
elseif ~strcmpi(x.M0,'UseControlAsM0') && isempty(regexpi(x.Q.Sequence, 'spiral'))
    % only register if the M0 and mean control are not identical
    % Which they are when there is no separate M0, but ASL was
    % acquired without background suppression & the mean control image
    % is used as M0 image

    % also skip registration of the M0 for 3D spiral, too poor resolution
    % and might worsen registration. PM: This part could be improved for 3D spiral

    if  xASL_exist(x.P.Path_mean_control,'file')
        refPath       = x.P.Path_mean_control;
    else % register M0 to T1w
        PathMNIMask   = fullfile(x.D.MapsSPMmodifiedDir, 'brainmask.nii'); % MNI brainmask
        xASL_im_SkullStrip(x.P.Path_T1, PathMNIMask, x);

        refPath       = x.P.Path_mask_T1;
        xASL_spm_smooth(refPath, [5 5 5],refPath);
    end



    %% 1A) Clip, mask & equalize contrast & intensity range
    %     xASL_im_EqualizeContrastImages( x.P.Path_rrM0, refIM );

    % First do the center of mass alignment
    if x.settings.bAutoACPC
        xASL_im_CenterOfMass(x.P.Path_rM0, {x.P.Path_M0} );
    end

    %% 2B) remove biasfields
    xASL_spm_BiasfieldCorrection(x.P.Path_rM0, x.D.SPMDIR, x.settings.Quality, [], x.P.Path_rrM0);
    xASL_spm_BiasfieldCorrection(refPath, x.D.SPMDIR, x.settings.Quality, [], refPath);


    %% 3C) Rigid-body registration
    xASL_spm_coreg(refPath, x.P.Path_rrM0, {x.P.Path_M0;x.P.Path_rM0}, x);
    xASL_delete(x.P.Path_mask_T1);
    xASL_delete(x.P.Path_rrM0);
end


%% -----------------------------------------------------------------------------------------------
%% 2) Quantify M0 (correction scale slope & incomplete T1 recovery)
M0_im = xASL_quant_M0(x.P.Path_rM0, x);



%% -----------------------------------------------------------------------------------------------
%% 3A) Conventional M0 masking & minor smoothing (doesnt work with smooth ASL images)
if x.settings.M0_conventionalProcessing
    % Conventional M0 processing, should be performed in native space
    % We
    % 1) perform the processing & masking in native space
    % 2) interpolate to standard space & mask again in standard space
    % Since this conventional method works with voxel-wise contrast differences, it needs to be
    % performed in native space. The transformation to standard space may play with the masking,
    % hence why we do the masking in both spaces

    fprintf('%s\n','Running conventional M0 processing method');

    xASL_spm_reslice( x.P.Path_ASL4D, x.P.Path_rM0, [], [], x.settings.Quality, x.P.Path_rM0, 1 ); % make sure M0 is in ASL space
    M0_nii          = xASL_io_ReadNifti( x.P.Path_rM0);
    x.VoxelSize     = M0_nii.hdr.pixdim(2:4);
    M0_im           = xASL_im_ProcessM0Conventional(M0_im, x); % also masks in native space
    % save image & mask
    xASL_io_SaveNifti(x.P.Path_rM0, x.P.Path_rM0, M0_im,[],0);
    xASL_io_SaveNifti(x.P.Path_rM0, x.P.Path_mask_M0, M0_im>0,8,0);

    % also transform to standard space
    InList  = {x.P.Path_rM0;x.P.Path_mask_M0};
    OutList = {x.P.Pop_Path_M0;x.P.Pop_Path_mask_M0};
    
    if exist(x.P.Path_mean_PWI_Clipped_sn_mat, 'file') % Backwards compatability, and also needed for the Affine+DCT co-registration of ASL-T1w
        AffineTransfPath = x.P.Path_mean_PWI_Clipped_sn_mat;
    else
        AffineTransfPath = [];
    end

    xASL_spm_deformations(x, InList, OutList, 1, [], AffineTransfPath, x.P.Path_y_ASL);

    % mask in standard space (native space masking already done)
	maskTmp = xASL_io_Nifti2Im(x.P.Pop_Path_M0);
    maskIM  = maskTmp.* (xASL_io_Nifti2Im(x.P.Pop_Path_mask_M0)==1);
	maskIM(maskTmp == 0) = NaN;
    xASL_io_SaveNifti(x.P.Pop_Path_M0,x.P.Pop_Path_M0,maskIM,8,0);

    % delete temporary files
    xASL_delete(x.P.Path_mask_M0);
    xASL_delete(x.P.Pop_Path_mask_M0);

    % Copy for visualization
    xASL_Copy(x.P.Pop_Path_M0,x.P.Pop_Path_noSmooth_M0,1);

%% -----------------------------------------------------------------------------------------------
%% 3B) New ExploreASL strategy for M0 masking & smoothing
%   Currently, this is performed in standard space, but this can be
%   optimized by running this in native space
else
    % Save image
    xASL_io_SaveNifti(x.P.Path_rM0,x.P.Path_rM0,M0_im,[],0);
    % Transform to standard space

    if exist(x.P.Path_mean_PWI_Clipped_sn_mat, 'file') % Backwards compatability, and also needed for the Affine+DCT co-registration of ASL-T1w
        AffineTransfPath = x.P.Path_mean_PWI_Clipped_sn_mat;
    else
        AffineTransfPath = [];
    end

    xASL_spm_deformations(x, x.P.Path_rM0, x.P.Pop_Path_M0, 1, [], AffineTransfPath, x.P.Path_y_ASL);
    % Copy for visualization (before smoothing/masking)
    xASL_Copy(x.P.Pop_Path_M0,x.P.Pop_Path_noSmooth_M0,1);

    fprintf('%s\n','Running new M0 processing method');
    % run the new M0 image processing
    IM = xASL_im_M0ErodeSmoothExtrapolate(xASL_io_Nifti2Im(x.P.Pop_Path_M0), x);
    xASL_io_SaveNifti(x.P.Pop_Path_M0, x.P.Pop_Path_M0, IM);
    % Copy M0 biasfield to native space for native space quantification

    if exist(x.P.Path_mean_PWI_Clipped_sn_mat, 'file') % Backwards compatability, and also needed for the Affine+DCT co-registration of ASL-T1w
        AffineTransfPath = x.P.Path_mean_PWI_Clipped_sn_mat;
    else
        AffineTransfPath = [];
    end

    xASL_spm_deformations(x, x.P.Pop_Path_M0, x.P.Path_rM0, 1, x.P.Path_PWI, AffineTransfPath, x.P.Path_y_ASL);
    % this last step also ensures that x.P.Path_rM0 is resliced to x.P.Path_PWI
end

%%
xASL_adm_CheckFileCount(x.D.PopDir, [x.P.M0 '_' x.P.SubjectID '_' x.P.SessionID '.nii'], 1);


end
