function xASL_wrp_RegisterASL(x)
%xASL_wrp_RegisterASL Submodule of ExploreASL ASL Module, that registers
%ASL to T1w (or potentially other structural images)
%
% FORMAT: xASL_wrp_RegisterASL(x)
%
% INPUT:
%   x  - structure containing fields with all information required to run this submodule (REQUIRED)
%
% OUTPUT: n/a (registration changes the NIfTI orientation header only,
%              with the exception of the affine transformation, which is
%              saved separately as x.P.Path_mean_PWI_Clipped_sn_mat
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This submodule registers ASL images to T1w space, by using a
% combination of the following registration techniques:
%
% M0-T1w rigid-body  -> this works well in 2D EPI sequences
% PWI-pGM rigid-body -> this is robust across sequences with different
%                       readouts and consequently different effective spatial resolutions. With
%                       low spatial resolution (e.g. GE 3D spiral product sequence), M0-T1w
%                       registration may not work, but PWI-pGM will work.
%                       PWI-pGM registration fails with large (vascular) artifacts, therefore
%                       this is performed only with relatively low spatial CoV.
% PWI-pGM affine     -> If the spatial CoV is sufficiently low, this can
%                       improve the registration
%
% These images are registered to ASL templates that were inversely
% transformed from MNI to the T1w space (& resampled to the ASL space)
% As this would have an ever higher similarity with the M0 & PWI
%
% This submodule performs the following steps:
%
% 1)    Create template images in native space
% 2)    Optimizes ASL image contrast for registration
% 3)    Registration ASL -> anat
%         x.bAffineRegistration - boolean to specify if the ASL-T1w rigid-body
%                                 registration is followed up by an affine
%                                 registration (OPTIONAL, DEFAULT = 2)
%                          - 0 = affine registration disabled
%                          - 1 = affine registration enabled
%                          - 2 = affine registration automatically chosen based on
%                                spatial CoV of PWI
%
% EXAMPLE: xASL_wrp_RegisterASL(x);
% __________________________________
% Copyright (C) 2015-2019 ExploreASL



%% ----------------------------------------------------------------------------------------
%% 0)   Administration

%% Use either original or motion estimated ASL4D
% Use despiked ASL only if spikes were detected and new file has been created
% Otherwise, despiked_raw_asl = same as original file
if ~xASL_exist(x.P.Path_despiked_ASL4D,'file')
    x.P.Path_despiked_ASL4D = x.P.Path_ASL4D;
end

% Remove pre-existing affine transformations, if we repeat registration
xASL_delete(x.P.Path_mean_PWI_Clipped_sn_mat);

%%      First run CenterOfMass registration & create ASL-standard space flow field
%       If no T1 flow field exists, create an identity flowfield.
%       We now assume the structural module hasn't run, and we simply want to run the ASL module to quickly check how the images look like
%       So we only run the automatic Center of Mass ACPC alignment
if x.bAutoACPC
    xASL_im_CenterOfMass(x.P.Path_despiked_ASL4D,{x.P.Path_M0});
end

if ~xASL_exist(x.P.Path_y_T1,'file')
    warning('Didnt find a structural scan, using MNI registration instead!!!!!!!!!!!' );
    IDmatrixPath    = fullfile(x.D.MapsSPMmodifiedDir,'Identity_Deformation_y_T1.nii');
    xASL_Copy(IDmatrixPath, x.P.Path_y_ASL,1);
else
    xASL_wrp_CreateASLDeformationField(x, true); % Smooth T1 deformation field into ASL resolution
end

%% ----------------------------------------------------------------------------------------
%% 1)   Create template images in native space

% Same for all sequences
Bias_Native         = fullfile(x.SESSIONDIR,'ATT_BiasField.nii');
Bias_MNI            = fullfile(x.D.TemplateDir,'ATT_BiasField.nii');
Vasc_MNI            = fullfile(x.D.TemplateDir,'MaxVesselTemplate.nii');
Vasc_Native         = fullfile(x.SESSIONDIR,'VascularArtifact_Template.nii');

Mean_Native         = fullfile(x.SESSIONDIR,'Mean_CBF_Template.nii');
Mask_Native         = fullfile(x.SESSIONDIR,'Mask_Template.nii');
raw_Native          = fullfile(x.SESSIONDIR,'RawTemplate.nii');

% Differs between sequences
if      strcmp(x.Sequence,'2D_EPI') && ~isempty(regexp(x.Vendor,'Philips'))
        Mean_MNI            = fullfile(x.D.TemplateDir,'Philips_2DEPI_Bsup_CBF.nii');
        Mask_MNI            = fullfile(x.D.TemplateDir,'Philips_2DEPI_Bsup_QC_mask.nii');

        if  x.Q.BackGrSupprPulses==0
            % No background suppression
            raw_MNI             = fullfile(x.D.TemplateDir,'Philips_2DEPI_noBsup_Control.nii');
        else % background suppression
            raw_MNI             = fullfile(x.D.TemplateDir,'Philips_2DEPI_Bsup_CBF.nii');
        end

elseif  strcmp(x.Sequence,'2D_EPI') && ~isempty(regexp(x.Vendor,'Siemens'))
        Mean_MNI            = fullfile(x.D.TemplateDir,'Siemens_2DEPI_PCASL_noBsup_CBF.nii');
        Mask_MNI            = fullfile(x.D.TemplateDir,'Philips_2DEPI_Bsup_QC_mask.nii');
        raw_MNI             = fullfile(x.D.TemplateDir,'Siemens_2DEPI_PCASL_noBsup_Control.nii');

        % in 3D readouts, background suppression is on by default, but this
        % doesn't matter for the average unsubtracted image, because there
        % is no slice-wise background suppression efficiency gradient
        % suspected

elseif  strcmp(x.Sequence,'3D_GRASE')
        raw_MNI             = fullfile(x.D.TemplateDir,'Siemens_3DGRASE_PCASL_Control_BiasfieldCorr_MoodStudy.nii');
        Mean_MNI            = fullfile(x.D.TemplateDir,'Siemens_3DGRASE_PASL_CBF.nii');
        Mask_MNI            = fullfile(x.D.TemplateDir,'Siemens_3DGRASE_PASL_QC_mask.nii');
elseif  strcmp(x.Sequence,'3D_spiral')
        Mean_MNI            = fullfile(x.D.TemplateDir,'GE_3Dspiral_Product_CBF.nii');
        Mask_MNI            = fullfile(x.D.MapsSPMmodifiedDir,'ParenchymNarrow.nii');
else
        error('Unknown sequence/readout');
end

if (~exist(Mean_Native,'file') || ~exist(Mask_Native,'file') || ~exist(raw_Native,'file') || ~exist(Vasc_Native,'file')) || strcmp(x.P.SessionID,x.SESSIONS{1})
    % We assume the same sequence for all sessions, so only have to warp
    % these templates to native space once. So we check if they exist, or
    % if we are in the first session (then we redo this)

    % For GE 3D spiral, skip raw image
    if strcmp(x.Sequence,'3D_spiral')
        xASL_spm_deformations(x,{Mean_MNI;Bias_MNI;Mask_MNI},{Mean_Native;Bias_Native;Mask_Native},1,x.P.Path_despiked_ASL4D,[],x.P.Path_y_ASL);
    else
        % trilinear interpolation is fine for smooth template
        xASL_spm_deformations(x,{Mean_MNI;Bias_MNI;Vasc_MNI;Mask_MNI;raw_MNI},{Mean_Native;Bias_Native;Vasc_Native;Mask_Native;raw_Native},1,x.P.Path_despiked_ASL4D,[],x.P.Path_y_ASL);
    end
end




%% ----------------------------------------------------------------------------------------
%% 2)    Optimizes ASL image contrast for registration
% Here we create a temporary dummy ASL image of which the image contrast is curated,
% for registration only
xASL_io_PairwiseSubtraction(x.P.Path_despiked_ASL4D, x.P.Path_mean_PWI_Clipped, 0, 0, x); % create PWI & mean_control
% Clip & compress the image, deal with contrast used for registration
tIM                         = xASL_im_ClipExtremes(x.P.Path_mean_PWI_Clipped, 0.95, 0.6); % careful, this cannot be rerun, once
tIM(tIM==min(tIM(:)))       = 0; % minimal intensities are set to 0
xASL_io_SaveNifti(x.P.Path_mean_PWI_Clipped, x.P.Path_mean_PWI_Clipped, tIM, [], 0);




%% ----------------------------------------------------------------------------------------
%% 3)    Registration ASL -> anat
% Here a temporary CBF image is created, which will be used for
% registration to T1 GM prob map & can be used for registration to previous ASL sessions.
% PWI-based registration can be preferable over EPI-based registration if
% background suppression and/or other 3D readout techniques are used.

% Here we first create a mask
% First check the initial alignment, otherwise first register with template
xASL_spm_reslice(Mask_Native, x.P.Path_mean_PWI_Clipped, x.P.Path_mean_PWI_Clipped_sn_mat, 0, x.Quality, x.P.Path_rmean_PWI_Clipped, 1 );
MaskASL     = xASL_im_ConvertMap2Mask(xASL_io_Nifti2Im(x.P.Path_rmean_PWI_Clipped));
MaskTemplate= logical(xASL_io_Nifti2Im(Mask_Native));
xASL_delete(x.P.Path_rmean_PWI_Clipped);

% compute dice coefficient
imA         = MaskASL;
imB         = MaskTemplate;
DiceCoeff   = xASL_im_ComputeDice(imA,imB);


% DiceCoeff should be around 0.5-0.6 after all registration, depending on
% correct mask/image comparison

OtherList = '';
OtherList{1,1} = x.P.Path_despiked_ASL4D;
OtherList{end+1,1} = x.P.Path_PWI;
OtherList{end+1,1} = x.P.Path_M0;
OtherList{end+1,1} = x.P.Path_ASL4D_RevPE;

if DiceCoeff(1)<0.4 && x.Quality && xASL_exist(raw_Native,'file')
    fprintf('%s\n','Poor initial registration, first registering with template');

    OtherList{end+1,1} = x.P.Path_mean_control;

    xASL_spm_coreg(raw_Native, x.P.Path_mean_PWI_Clipped, OtherList, x, [6 3]);
end

if xASL_exist(x.P.Path_mean_control,'file')

    OtherList{4,1} = x.P.Path_mean_PWI_Clipped;

    if x.Quality
        xASL_spm_coreg(raw_Native, x.P.Path_mean_control, OtherList, x);
    else
        xASL_spm_coreg(raw_Native, x.P.Path_mean_control, OtherList, x,2);  % low quality registration
    end

    spatCoVit = xASL_im_GetSpatialCovNativePWI(x);
    if  spatCoVit>0.6 && x.bPWIRegistration
        nIT                         = 0;
        fprintf('%s\n','High spatial CoV, skipping CBF-based registration');
    else
        nIT                         = 1; % perform 1 PWI registration, first was with raw image
    end

else
    spatCoVit(1)                    = 0.5; % safe start
    nIT                             = 2; % perform 2 PWI registrations
end

if ~x.Quality && nIT>1
    nIT = 1; % speed up for low quality
end

if nIT>0
    xASL_im_CreatePseudoCBF(x,spatCoVit(1));
    % keep this same for all sequences, 3D spiral will simply have a lower spatial CoV because of smoothness

    for iT=1:nIT
        OtherList{4,1}                  = x.P.Path_mean_control; % rest of otherlist has been defined above
        xASL_spm_coreg(x.P.Path_PseudoCBF, x.P.Path_mean_PWI_Clipped, OtherList, x);

        spatCoV                         = xASL_im_GetSpatialCovNativePWI(x);
        xASL_im_CreatePseudoCBF(x,spatCoV);
        spatCoVit(iT+1)                 = spatCoV;
    end

    %% Affine registration
    if isfield(x,'bAffineRegistration') && x.bAffineRegistration<2
        bAffineRegistration = x.bAffineRegistration;
    else % if not set, default here is to only do affine regitration for high quality processing & low spatial CoV
        bAffineRegistration = x.Quality && spatCoVit(end)<0.4;
    end

    if bAffineRegistration % perform affine registration
        xASL_spm_affine(x.P.Path_mean_PWI_Clipped, x.P.Path_PseudoCBF, 5, 5, OtherList);
    end
end


%% ----------------------------------------------------------------------------------------
%% Delete temporary files
if x.DELETETEMP
    File2Del  = {Mean_Native Bias_Native Vasc_Native Mask_Native raw_Native x.P.Path_mean_PWI_Clipped x.P.Path_mean_control};
    for iL=1:length(File2Del)
        xASL_delete(File2Del{iL});
    end
end

fprintf('\n%s\n','--------------------------------------------------------------------');

fprintf('%s\n',[num2str(length(spatCoVit)) ' registration iterations:']);
for iT=1:length(spatCoVit)
    fprintf('%s\n',['Iteration ' num2str(iT) ', spatial CoV = ' num2str(100*spatCoVit(iT),3) '%']);
end
fprintf('%s\n\n','--------------------------------------------------------------------');


end
