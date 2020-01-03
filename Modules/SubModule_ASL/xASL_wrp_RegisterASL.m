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
% 0)    Administration: 
%     - A) ASL4D is dealth with, if motion peaks were removed this is called
%       despiked_ASL4D
%     - B) a default "OtherList" is specified. This is used every
%         registration instance, except for removing the ref and src NIfTIs
%         used in the registration instance. Also, inside the registration
%         function the unexisting OtherList NIfTIs are skipped
%     - C) Define paths to the ASL templates
%     - D) Previous registration output files are removed
%     - E) native->MNI transformation flow field y_T1.nii is smoothed to the 
%          effective ASL resolution y_ASL.nii
%     - F) Registration contrasts are dealth with:
%       x.bRegistrationContrast - specifies the image contrast used for
%                                 registration (OPTIONAL, DEFAULT = 2):
%                           - 0 = Control-T1w
%                           - 1 = CBF - pseudoCBF from template/pGM+pWM
%                           - 2 = automatic (mix of both)
%     - G) Dummy src NIfTIs are created:
%       mean_control.nii to register with T1w
%       mean_PWI_Clipped.nii to register with pseudoCBF
%
% 1)    Registration Center of Mass
% 2)    Registration ASL -> anat (Control->T1w)
% 3)    Registration CBF->pseudoCBF


%       x.bAffineRegistration - specifies the ASL-T1w rigid-body
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
% A) Use either original or motion estimated ASL4D
% Use despiked ASL only if spikes were detected and new file has been created
% Otherwise, despiked_raw_asl = same as original file
if ~xASL_exist(x.P.Path_despiked_ASL4D,'file')
    x.P.Path_despiked_ASL4D = x.P.Path_ASL4D;
end

%% B) Manage OtherList
% Define OtherList for registration
% Here we mention all possible files that need to be in registration
% All functions below will remove those that are unexisting, or used in the
% registration estimation.
BaseOtherList = {x.P.Path_despiked_ASL4D x.P.Path_mean_control x.P.Path_M0 x.P.Path_PWI x.P.Path_mean_PWI_Clipped x.P.Path_ASL4D_RevPE};

if ~strcmp(x.P.Path_despiked_ASL4D, x.P.Path_ASL4D)
    BaseOtherList{end+1} = x.P.Path_ASL4D; % keep original ASL4D aligned as well
end

%% C) Define paths to the ASL templates
% Same for all sequences
x.Bias_Native = fullfile(x.SESSIONDIR,'ATT_BiasField.nii');
x.Bias_MNI = fullfile(x.D.TemplateDir,'ATT_BiasField.nii');
x.Vasc_Native = fullfile(x.SESSIONDIR,'VascularArtifact_Template.nii');
x.Vasc_MNI = fullfile(x.D.TemplateDir,'MaxVesselTemplate.nii');

x.Mean_Native = fullfile(x.SESSIONDIR,'Mean_CBF_Template.nii');
x.Mask_Native = fullfile(x.SESSIONDIR,'Mask_Template.nii');
x.raw_Native = fullfile(x.SESSIONDIR,'RawTemplate.nii');

x.PathMask = fullfile(x.SESSIONDIR, 'MaskASL.nii');
x.Path_PseudoTissue = fullfile(x.SESSIONDIR, 'PseudoTissue.nii');

% Differs between sequences
if      strcmp(x.Sequence,'2D_EPI') && ~isempty(regexp(x.Vendor,'Philips'))
        x.Mean_MNI = fullfile(x.D.TemplateDir,'Philips_2DEPI_Bsup_CBF.nii');
        x.Mask_MNI = fullfile(x.D.TemplateDir,'Philips_2DEPI_Bsup_QC_mask.nii');
        x.raw_MNI = fullfile(x.D.TemplateDir,'Philips_2DEPI_noBsup_Control.nii');

elseif  strcmp(x.Sequence,'2D_EPI') && ~isempty(regexp(x.Vendor,'Siemens'))
        x.Mean_MNI = fullfile(x.D.TemplateDir,'Siemens_2DEPI_PCASL_noBsup_CBF.nii');
        x.Mask_MNI = fullfile(x.D.TemplateDir,'Philips_2DEPI_Bsup_QC_mask.nii');
        x.raw_MNI = fullfile(x.D.TemplateDir,'Siemens_2DEPI_PCASL_noBsup_Control.nii');

        % in 3D readouts, background suppression is on by default, but this
        % doesn't matter for the average unsubtracted image, because there
        % is no slice-wise background suppression efficiency gradient
        % suspected

        % The raw_MNI template was used to register the Control to,
        % but Control->T1w seems to outperform this Control template (it is
        % too smooth). The raw_MNI was always skipped for 3D_spiral, as
        % this GE sequence usually is provided without separate
        % control-label images
        
elseif  strcmp(x.Sequence,'3D_GRASE')
        x.raw_MNI = fullfile(x.D.TemplateDir,'Siemens_3DGRASE_PCASL_Control_BiasfieldCorr_MoodStudy.nii');
        x.Mean_MNI = fullfile(x.D.TemplateDir,'Siemens_3DGRASE_PASL_CBF.nii');
        x.Mask_MNI = fullfile(x.D.TemplateDir,'Siemens_3DGRASE_PASL_QC_mask.nii');
elseif  strcmp(x.Sequence,'3D_spiral')
        x.Mean_MNI = fullfile(x.D.TemplateDir,'GE_3Dspiral_Product_CBF.nii');
        x.Mask_MNI = fullfile(x.D.MapsSPMmodifiedDir,'ParenchymNarrow.nii');
else
        error('Unknown sequence/readout');
end


%% D) Remove pre-existing registration information, if we repeat registration
xASL_delete(x.P.Path_mean_PWI_Clipped_sn_mat);
xASL_delete(x.P.Path_mean_PWI_Clipped);
xASL_delete(x.P.Path_mean_control);
if strcmp(x.P.SessionID,'ASL_1') || x.nSessions==1
    xASL_delete(x.Path_PseudoTissue);
    xASL_delete(x.Bias_Native);
    xASL_delete(x.Vasc_Native);
    xASL_delete(x.Mask_Native);
    xASL_delete(x.Mean_Native);
    xASL_delete(x.raw_Native);
    xASL_delete(x.P.Path_PseudoCBF);
    xASL_delete(fullfile(x.SESSIONDIR,'MaskASL.nii'));
end
 
    
    
%% E) % Smooth T1 deformation field into ASL resolution
% If no T1 flow field exists, create an identity flowfield
% So we can still process ASL images without the
% structural module
if ~xASL_exist(x.P.Path_y_ASL,'file') || strcmp(x.P.SessionID,'ASL_1') || x.nSessions==1
    if ~xASL_exist(x.P.Path_y_T1,'file')
        warning('Didnt find a structural scan, using MNI registration instead!!!!!!!!!!!');
        IDmatrixPath = fullfile(x.D.MapsSPMmodifiedDir,'Identity_Deformation_y_T1.nii');
        xASL_Copy(IDmatrixPath, x.P.Path_y_ASL, true);
    else
        xASL_wrp_CreateASLDeformationField(x, true);
    end
end

%% F) Manage registration contrasts that we will use
if isfield(x,'bRegistrationContrast') && x.bRegistrationContrast==0
    bRegistrationControl = true;
    bRegistrationCBF = false;
elseif isfield(x,'bRegistrationContrast') && x.bRegistrationContrast==1
    bRegistrationControl = false;
    bRegistrationCBF = true;
else
    bRegistrationControl = true;
    bRegistrationCBF = true;
end

%% G) Here we create a temporary dummy ASL image of which the image contrast is curated,
% for registration only
xASL_io_PairwiseSubtraction(x.P.Path_despiked_ASL4D, x.P.Path_mean_PWI_Clipped, 0, 0, x); % create PWI & mean_control
if bRegistrationCBF
    % Clip & compress the image, deal with contrast used for registration
    tIM = xASL_im_ClipExtremes(x.P.Path_mean_PWI_Clipped, 0.95, 0.6); % careful, this cannot be rerun, once
    tIM(tIM==min(tIM(:))) = 0; % minimal intensities are set to 0
    xASL_io_SaveNifti(x.P.Path_mean_PWI_Clipped, x.P.Path_mean_PWI_Clipped, tIM, [], 0);
end



%% ----------------------------------------------------------------------------------------
%% 1)   Registration CenterOfMass
% We now assume the structural module hasn't run, and we simply want to run the ASL module to quickly check how the images look like
% So we only run the automatic Center of Mass ACPC alignment
if x.bAutoACPC
    OtherList = xASL_adm_RemoveFromOtherList(BaseOtherList, {x.P.Path_despiked_ASL4D});
    % x.P.Path_despiked_ASL4D is padded to the end of the list
    xASL_im_CenterOfMass(x.P.Path_despiked_ASL4D,OtherList);
end

    
%% ----------------------------------------------------------------------------------------
%% 2)    Registration Control->T1w
% Here a temporary CBF image is created, which will be used for
% registration to T1 GM prob map & can be used for registration to previous ASL sessions.
% PWI-based registration can be preferable over EPI-based registration if
% background suppression and/or other 3D readout techniques are used.

% Here we first create a mask
% First check the initial alignment, otherwise first register with template
% xASL_spm_reslice(Mask_Native, x.P.Path_mean_PWI_Clipped, x.P.Path_mean_PWI_Clipped_sn_mat, 0, x.Quality, x.P.Path_rmean_PWI_Clipped, 1);
% MaskASL = xASL_im_ConvertMap2Mask(xASL_io_Nifti2Im(x.P.Path_rmean_PWI_Clipped));
% MaskTemplate= logical(xASL_io_Nifti2Im(Mask_Native));
% xASL_delete(x.P.Path_rmean_PWI_Clipped);

% % compute dice coefficient:
% % DiceCoeff should be around 0.5-0.6 after all registration, depending on
% % correct mask/image comparison
% imA = MaskASL;
% imB = MaskTemplate;
% DiceCoeff = xASL_im_ComputeDice(imA,imB);

if bRegistrationControl
% 1) Initial mean control registrations, if available
    if ~xASL_exist(x.P.Path_mean_control,'file')
        warning('Requested control-T1w registration but Mean_control.nii didnt exist');
    else
        OtherList = xASL_adm_RemoveFromOtherList(BaseOtherList, {x.P.Path_mean_control});
        if ~x.Quality
            xASL_spm_coreg(x.P.Path_T1, x.P.Path_mean_control, OtherList, x,[9 6]);
        else
            xASL_spm_coreg(x.P.Path_T1, x.P.Path_mean_control, OtherList, x);
        end
    end
end



%% ----------------------------------------------------------------------------------------
%% 3)   Registration CBF->pseudoCBF
if bRegistrationCBF
    % This creates the reference images, the downsampled pseudoTissue
    % (pGM+pWM) as well as the native space copies of templates for CBF,
    % ATT biasfield and vascular peaks
    xASL_im_CreatePseudoCBF(x, 0);


    spatCoVit = xASL_im_GetSpatialCovNativePWI(x);
    if spatCoVit>0.6
        nIT = 0;
        fprintf('%s\n','High spatial CoV, skipping CBF-based registration');
    elseif ~x.Quality
        nIT = 1; % speed up for low quality
    elseif bRegistrationControl
        nIT = 1; % if control-T1w registration is already performed, 1 iteration is OK here
    else
        nIT = 2;
    end


    % 2) Repeat CBF registrations, with iteratively better estimate of the
    % vascular/tissue perfusion ratio of the template
    if nIT>0
        % keep this same for all sequences, 3D spiral will simply have a lower spatial CoV because of smoothness
        for iT=1:nIT
            xASL_im_CreatePseudoCBF(x, spatCoVit(1));
            OtherList = xASL_adm_RemoveFromOtherList(BaseOtherList, {x.P.Path_mean_PWI_Clipped});
            xASL_spm_coreg(x.P.Path_PseudoCBF, x.P.Path_mean_PWI_Clipped, OtherList, x);

            spatCoVit(iT+1) = xASL_im_GetSpatialCovNativePWI(x);
        end

        %% Affine registration
        if isfield(x,'bAffineRegistration') && x.bAffineRegistration<2
            bAffineRegistration = x.bAffineRegistration;
        else % if not set, default here is to only do affine regitration for high quality processing & low spatial CoV
            bAffineRegistration = x.Quality && spatCoVit(end)<0.4;
        end

        if bAffineRegistration % perform affine registration
            xASL_im_CreatePseudoCBF(x, spatCoVit(end));
            % apply also to mean_PWI_clipped and other files
            xASL_spm_affine(x.P.Path_mean_PWI_Clipped, x.P.Path_PseudoCBF, 5, 5, BaseOtherList);
            spatCoVit(iT+2) = xASL_im_GetSpatialCovNativePWI(x);
        end
    end
    
    fprintf('\n%s\n','--------------------------------------------------------------------');

    fprintf('%s\n',[num2str(length(spatCoVit)) ' registration iterations:']);
    for iT=1:length(spatCoVit)
        fprintf('%s\n',['Iteration ' num2str(iT) ', spatial CoV = ' num2str(100*spatCoVit(iT),3) '%']);
    end
    fprintf('%s\n\n','--------------------------------------------------------------------');
end


%% ----------------------------------------------------------------------------------------
%% Delete temporary files
if x.DELETETEMP
    File2Del = {x.Mean_Native x.Bias_Native x.Vasc_Native x.Mask_Native x.raw_Native x.P.Path_mean_PWI_Clipped x.P.Path_mean_control x.P.Path_PseudoCBF x.PathMask};
    for iL=1:length(File2Del)
        xASL_delete(File2Del{iL});
    end
end


end


function [OtherList] = xASL_adm_RemoveFromOtherList(BaseOtherList, List2Remove)
%xASL_adm_RemoveFromOtherList Take the base "OtherList" and remove the
%files that are provided as second input argument, in cell format, from the
%BaseOtherList, outputting as OtherList that can be used by any image
%processing, to apply estimations/transformations to. Note that inside
%these functions, most likely the unexisting files will be removed from the
%OtherList as well

OtherList = BaseOtherList;

for iList=1:length(List2Remove)
    IndexIs = find(cellfun(@(y) strcmp(y,List2Remove{iList}), OtherList));
    if IndexIs==1
        IndexList = [2:length(OtherList)];
    elseif IndexIs==length(OtherList)
        IndexList = [1:length(OtherList)-1];
    else
        IndexList = [1:IndexIs-1 IndexIs+1:length(OtherList)];
    end
    OtherList = OtherList(IndexList);
end

end