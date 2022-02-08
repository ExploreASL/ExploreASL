function x = xASL_wrp_RegisterASL(x)
%xASL_wrp_RegisterASL Submodule of ExploreASL ASL Module, that registers
%ASL to T1w (or potentially other structural images)
%
% FORMAT: xASL_wrp_RegisterASL(x)
%
% INPUT:
%   x  - structure containing fields with all information required to run this submodule (REQUIRED)
%
% OUTPUT: x - updated structure (mainly the Tanimoto coefficient of the final registration),
%             and the registration also changes the NIfTI orientation header
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DEVELOPER:
% PM: the vascular template registration may need some improvement
% PM: this function can be divided into subfunctions for readability and to be less bug-prone
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This submodule registers ASL images to T1w space, by using a
% combination of the registration techniques below. Note that in the
% absence of raw structural files (i.e. T1.nii[.gz] or T1_ORI.nii[.gz],
% it will recreate dummy files from standard space to do this registration
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
% 0.    Administration:
%     - A. ASL4D is dealth with, if motion peaks were removed this is called
%       despiked_ASL4D
%     - B. a default "OtherList" is specified. This is used every
%         registration instance, except for removing the ref and src NIfTIs
%         used in the registration instance. Also, inside the registration
%         function the unexisting OtherList NIfTIs are skipped
%     - C. Define paths to the ASL templates
%     - D. Previous registration output files are removed
%     - E. Allow registration without structural data
%     - F. native->MNI transformation flow field y_T1.nii is smoothed to the
%          effective ASL resolution y_ASL.nii
%     - G. Registration contrasts are dealth with:
%       x.modules.asl.bRegistrationContrast - specifies the image contrast used for
%                                 registration (OPTIONAL, DEFAULT = 2):
%                           - 0 = Control->T1w
%                           - 1 = CBF->pseudoCBF from template/pGM+pWM
%                                 (skip if sCoV>0.667)
%                           - 2 = automatic (mix of both)
%                           - 3 = option 2 & force CBF->pseudoCBF irrespective of sCoV or Tanimoto coefficient
%     - H. Dummy src NIfTIs are created:
%          mean_control.nii to register with T1w
%          mean_PWI_Clipped.nii to register with pseudoCBF
%     - I. Create reference images, downsampled pseudoTissue
%
% 1.    Registration Center of Mass
% 2.    Registration ASL -> anat (Control->T1w)
%       (this step is only applied if it improves the Tanimoto coefficient)
% 3.    Registration CBF->pseudoCBF
%       (this step is only applied if it improves the Tanimoto coefficient). Also, this step is only
%       applied if the spatial CoV<0.67. Note that this is usually the case
%       for 3D scans because of their lower effective spatial resolution.
%
%       x.modules.asl.bAffineRegistration - specifies the ASL-T1w rigid-body registration is followed up by an affine
%                                 registration (OPTIONAL, DEFAULT = 0)
%                          - 0 = affine registration disabled
%                          - 1 = affine registration enabled
%                          - 2 = affine registration automatically chosen based on
%                                spatial CoV of PWI
%       x.modules.asl.bDCTRegistration - Specifies if to include the DCT registration on top of Affine, all other requirements for
%                            affine are thus also taken into account (OPTIONAL, DEFAULT = 0)
%                          - 0 = DCT registration disabled
%                          - 1 = DCT registration enabled if affine enabled and conditions for affine passed
%                          - 2 = DCT enabled as above, but use PVC on top of it to get the local intensity scaling right
%
% EXAMPLE: xASL_wrp_RegisterASL(x);
% __________________________________
% Copyright (C) 2015-2021 ExploreASL


%% ----------------------------------------------------------------------------------------
%% 0.   Administration
% A) Use either original or motion estimated ASL4D
% Use despiked ASL only if spikes were detected and new file has been created
% Otherwise, despiked_raw_asl = same as original file
if ~xASL_exist(x.P.Path_despiked_ASL4D,'file')
    x.P.Path_despiked_ASL4D = x.P.Path_ASL4D;
end

if ~isfield(x.modules.asl,'bRegistrationContrast') || isempty(x.modules.asl.bRegistrationContrast)
    x.modules.asl.bRegistrationContrast = 2; % register M0-T1w first, then CBF-pGM if sCoV<0.667
end

if ~isfield(x.modules.asl,'bAffineRegistration') || isempty(x.modules.asl.bAffineRegistration)
	x.modules.asl.bAffineRegistration = 0; % Default - affine disabled
end

% DCT off by default
if ~isfield(x.modules.asl,'bDCTRegistration') || isempty(x.modules.asl.bDCTRegistration)
	x.modules.asl.bDCTRegistration = 0;
end

% DCT runs only with affine, switching DCT on and affine off thus produces a warning
if x.modules.asl.bDCTRegistration && ~x.modules.asl.bAffineRegistration
	warning('DCT registration cannot run if affine is disabled');
end

% For DCT, force CBF<->pseudoCBF registration
if x.modules.asl.bDCTRegistration
	x.modules.asl.bRegistrationContrast = 3;
end

% By default, don't use dummy structural even if the structural image is missing
if ~isfield(x.modules.asl,'bUseMNIasDummyStructural') || isempty(x.modules.asl.bUseMNIasDummyStructural)
	x.modules.asl.bUseMNIasDummyStructural = false;
end
if ~isfield(x.modules.asl, 'RegistrationTemplate') || isempty(x.modules.asl.RegistrationTemplate)
     x.modules.asl.RegistrationTemplate = 'MNI_Structural'; % default
end


%% B. Manage OtherList
% Define OtherList for registration
% Here we mention all possible files that need to be in registration
% All functions below will remove those that are unexisting, or used in the
% registration estimation.
BaseOtherList = {x.P.Path_despiked_ASL4D, x.P.Path_mean_control, x.P.Path_M0, x.P.Path_PWI, x.P.Path_mean_PWI_Clipped, x.P.Path_ASL4D_RevPE,...
	x.P.Path_mean_PWI_Clipped_DCT,x.P.Path_mean_PWI_Clipped_ORI,...
    x.P.Path_ASL4D_ORI, fullfile(x.dir.SESSIONDIR, 'B0.nii'), fullfile(x.dir.SESSIONDIR, 'Unwarped.nii'), fullfile(x.dir.SESSIONDIR, 'Field.nii'), fullfile(x.dir.SESSIONDIR, 'TopUp_fieldcoef.nii')};

if ~strcmp(x.P.Path_despiked_ASL4D, x.P.Path_ASL4D)
    BaseOtherList{end+1} = x.P.Path_ASL4D; % keep original ASL4D aligned as well
end

%% C. Define paths to the ASL templates
% Same for all sequences
x.D.Bias_Native = fullfile(x.dir.SESSIONDIR,'ATT_BiasField.nii');
x.D.Bias_MNI = fullfile(x.D.TemplateDir,'ATT_BiasField.nii');
x.D.Vasc_Native = fullfile(x.dir.SESSIONDIR,'VascularArtifact_Template.nii');
x.D.Vasc_MNI = fullfile(x.D.TemplateDir,'MaxVesselTemplate.nii');

x.D.Mean_Native = fullfile(x.dir.SESSIONDIR,'Mean_CBF_Template.nii');
x.D.Mask_Native = fullfile(x.dir.SESSIONDIR,'Mask_Template.nii');
x.D.raw_Native = fullfile(x.dir.SESSIONDIR,'RawTemplate.nii');

x.D.PathMask = fullfile(x.dir.SESSIONDIR, 'MaskASL.nii');
x.D.Path_PseudoTissue = fullfile(x.dir.SESSIONDIR, 'PseudoTissue.nii');

% Differs between sequences
if x.bQASPERPhantom
    x.D.Mask_MNI = fullfile(x.D.TemplateDir,'QASPER', 'QASPER_QC_mask.nii');
    x.D.Mean_MNI = fullfile(x.D.TemplateDir,'Philips_2DEPI_Bsup_CBF.nii'); % needs fixing
    x.D.raw_MNI = fullfile(x.D.TemplateDir,'Philips_2DEPI_noBsup_Control.nii'); % needs fixing
    
elseif strcmpi(x.Q.Sequence,'2D_EPI') && ~isempty(regexpi(x.Q.Vendor,'Philips'))
        x.D.Mean_MNI = fullfile(x.D.TemplateDir,'Philips_2DEPI_Bsup_CBF.nii');
        x.D.Mask_MNI = fullfile(x.D.TemplateDir,'Philips_2DEPI_Bsup_QC_mask.nii');
        x.D.raw_MNI = fullfile(x.D.TemplateDir,'Philips_2DEPI_noBsup_Control.nii');

elseif  strcmpi(x.Q.Sequence,'2D_EPI') && ~isempty(regexpi(x.Q.Vendor,'(Siemens|GE)'))
        %% PM: quicky & dirty fix to run GE 2D EPI with the Siemens 2D EPI template
        % though the template choice may not have a significant effect, as
        % opposed to the inter-individual differences in geometric
        % distortion
        x.D.Mean_MNI = fullfile(x.D.TemplateDir,'Siemens_2DEPI_PCASL_noBsup_CBF.nii');
        x.D.Mask_MNI = fullfile(x.D.TemplateDir,'Philips_2DEPI_Bsup_QC_mask.nii');
        x.D.raw_MNI = fullfile(x.D.TemplateDir,'Siemens_2DEPI_PCASL_noBsup_Control.nii');

        % in 3D readouts, background suppression is on by default, but this
        % doesn't matter for the average unsubtracted image, because there
        % is no slice-wise background suppression efficiency gradient
        % suspected

        % The raw_MNI template was used to register the Control to,
        % but Control->T1w seems to outperform this Control template (it is
        % too smooth). The raw_MNI was always skipped for 3D_spiral, as
        % this GE sequence usually is provided without separate
        % control-label images

elseif  strcmpi(x.Q.Sequence,'3D_GRASE')
        x.D.raw_MNI = fullfile(x.D.TemplateDir,'Siemens_3DGRASE_PCASL_Control_BiasfieldCorr_MoodStudy.nii');
        x.D.Mean_MNI = fullfile(x.D.TemplateDir,'Siemens_3DGRASE_PASL_CBF.nii');
        x.D.Mask_MNI = fullfile(x.D.TemplateDir,'Siemens_3DGRASE_PASL_QC_mask.nii');
elseif  strcmpi(x.Q.Sequence,'3D_spiral')
        x.D.Mean_MNI = fullfile(x.D.TemplateDir,'GE_3Dspiral_Product_CBF.nii');
        x.D.Mask_MNI = fullfile(x.D.MapsSPMmodifiedDir,'ParenchymNarrow.nii');
else
        error('Unknown sequence/readout');
end


%% D. Remove pre-existing registration information, if we repeat registration
xASL_delete(x.P.Path_mean_PWI_Clipped_sn_mat);
xASL_delete(x.P.Path_mean_PWI_Clipped);
xASL_delete(x.P.Path_mean_control);
if strcmp(x.P.SessionID,'ASL_1') || x.dataset.nSessions==1
    xASL_delete(x.D.Path_PseudoTissue);
    xASL_delete(x.D.Bias_Native);
    xASL_delete(x.D.Vasc_Native);
    xASL_delete(x.D.Mask_Native);
    xASL_delete(x.D.Mean_Native);
    xASL_delete(x.D.raw_Native);
    xASL_delete(x.P.Path_PseudoCBF);
    xASL_delete(fullfile(x.dir.SESSIONDIR,'MaskASL.nii'));
end



%% E. Allow registration without structural data
StructuralDerivativesExist = xASL_exist(x.P.Path_y_T1, 'file') && xASL_exist(x.P.Path_c1T1, 'file') && xASL_exist(x.P.Path_c2T1, 'file');
StructuralRawExist = xASL_exist(x.P.Path_T1, 'file') || xASL_exist(x.P.Path_T1_ORI, 'file');

if StructuralRawExist && ~StructuralDerivativesExist
    error('Please run structural module first');
elseif ~StructuralRawExist && ~StructuralDerivativesExist
    if x.modules.asl.bUseMNIasDummyStructural
        xASL_wrp_UseTemplateAsDummyStructural(x, x.modules.asl.RegistrationTemplate);
    else
        error('Structural data missing, skipping ASL module; if this is undesired, set x.modules.asl.bUseMNIasDummyStructural=1');
    end
end

%% F. Smooth T1 deformation field into ASL resolution
% If no T1 flow field exists, create an identity flowfield
% So we can still process ASL images without the
% structural module
if ~xASL_exist(x.P.Path_y_ASL,'file') || strcmp(x.P.SessionID,'ASL_1') || x.dataset.nSessions==1
    xASL_im_CreateASLDeformationField(x, true);
end


%% G. Manage registration contrasts that we will use
if x.modules.asl.bContainsDeltaM
	bRegistrationControl = false;
    bRegistrationCBF = true;
elseif x.modules.asl.bRegistrationContrast==0
    bRegistrationControl = true;
    bRegistrationCBF = false;
elseif x.modules.asl.bRegistrationContrast==1
    bRegistrationControl = false;
    bRegistrationCBF = true;
else
    bRegistrationControl = true;
    bRegistrationCBF = true;
end

%% H. Here we create a temporary dummy ASL image of which the image contrast is curated,
% for registration only

% Load the image
ASL_im = xASL_io_Nifti2Im(x.P.Path_despiked_ASL4D);
if (size(ASL_im, 4) == 1) || x.modules.asl.bContainsDeltaM
	% Apparently, the subtraction was already done on the scanner/reconstruction
	
	% Save the mean of deltaM without subtracting with using the right coordinates
	xASL_io_SaveNifti(x.P.Path_despiked_ASL4D, x.P.Path_mean_PWI_Clipped, xASL_stat_MeanNan(ASL_im,4), 16, false);
elseif x.modules.asl.bTimeEncoded
	% Decoding of TimeEncoded data (Nifti is saved inside the function)
	ASL_im = xASL_quant_HadamardDecoding(x.P.Path_despiked_ASL4D, x.Q);
	
	% Hadamard Block size is calculated as number of TEs and the HadamardMatrixSize-1
	blockSize = x.Q.NumberEchoTimes * (x.Q.TimeEncodedMatrixSize-1);
	
	PWI = zeros(size(ASL_im,1), size(ASL_im,2), size(ASL_im,3), blockSize); % preallocate PWI
	
	for iBlock = 1:blockSize
		PWI(:,:,:,iBlock) = xASL_stat_MeanNan(ASL_im(:,:,:,iBlock:blockSize:end), 4); % Averaged PWI4D across repetitions
	end
	
	% Create single PWI for further steps in ASL module
	PWI = xASL_stat_MeanNan(PWI(:,:,:,1:x.Q.NumberEchoTimes:end),4); % Average across PLDs from each first TE
	
	% Save single PWI
	xASL_io_SaveNifti(x.P.Path_despiked_ASL4D, x.P.Path_mean_PWI_Clipped, PWI, 16, false);
else
	xASL_io_PairwiseSubtraction(x.P.Path_despiked_ASL4D, x.P.Path_mean_PWI_Clipped, 0, 0); % create PWI & mean_control
end

% MultiPLD: we can keep the PairwiseSubtraction function->because we have controls and labels for each PLD and take average of all PLDs
% (later on maybe we only need to use the later PLDs)

if bRegistrationCBF
    % This clipped_ORI is not used in the rest of the pipeline
	xASL_Copy(x.P.Path_mean_PWI_Clipped,x.P.Path_mean_PWI_Clipped_ORI,1);
    % Clip & compress the image, deal with contrast used for registration
    tIM = xASL_im_ClipExtremes(x.P.Path_mean_PWI_Clipped, 0.95, 0.6); % careful, this cannot be rerun, once
    tIM(tIM==min(tIM(:))) = 0; % minimal intensities are set to 0
    xASL_io_SaveNifti(x.P.Path_mean_PWI_Clipped, x.P.Path_mean_PWI_Clipped, tIM, [], 0);
end



%% I. Here we create the reference images, the downsampled pseudoTissue
% (pGM+pWM) as well as the native space copies of templates for CBF,
% ATT biasfield and vascular peaks
xASL_im_CreatePseudoCBF(x, 0);

TanimotoPerc = xASL_im_GetSpatialOverlapASL(x);


%% ----------------------------------------------------------------------------------------
%% 1.   Registration CenterOfMass
% We now assume the structural module hasn't run, and we simply want to run the ASL module to quickly check how the images look like
% So we only run the automatic Center of Mass ACPC alignment
if x.settings.bAutoACPC
    OtherList = xASL_adm_RemoveFromOtherList(BaseOtherList, {x.P.Path_despiked_ASL4D}); % x.P.Path_despiked_ASL4D is padded to the end of the list

    xASL_im_BackupAndRestoreAll(BaseOtherList, 1); % First backup all NIfTIs & .mat sidecars of BaseOtherList
    
    if x.bQASPERPhantom
        xASL_im_CenterOfMass(x.P.Path_despiked_ASL4D, OtherList, 0, [0;0;13.5]); 
        % QASPER needs different offset, it has a different center of mass than a human brain
    else
        xASL_im_CenterOfMass(x.P.Path_despiked_ASL4D, OtherList, 0); % default offset        
    end
    TanimotoPerc(end+1) = xASL_im_GetSpatialOverlapASL(x); % get new overlap score
    if TanimotoPerc(end)>=TanimotoPerc(end-1) % if alignment improved or remained same
        xASL_im_BackupAndRestoreAll(BaseOtherList, 3); % delete backup
    else % if alignment got worse
		TanimotoPerc = TanimotoPerc(1:end-1); % Remove the last Tanimoto number as this has been restored
        xASL_im_BackupAndRestoreAll(BaseOtherList, 2); % restore NIfTIs from backup
    end
end


%% ----------------------------------------------------------------------------------------
%% 2.    Registration Control->T1w
% Here we first create a mask
% First check the initial alignment, otherwise first register with template
% xASL_spm_reslice(Mask_Native, x.P.Path_mean_PWI_Clipped, x.P.Path_mean_PWI_Clipped_sn_mat, 0, x.settings.Quality, x.P.Path_rmean_PWI_Clipped, 1);
% MaskASL = xASL_im_ConvertMap2Mask(xASL_io_Nifti2Im(x.P.Path_rmean_PWI_Clipped));
% MaskTemplate= logical(xASL_io_Nifti2Im(Mask_Native));
% xASL_delete(x.P.Path_rmean_PWI_Clipped);

% compute Tanimoto coefficient:
% Tanimoto Coeff should be around 0.8 after all registration, depending on
% correct mask/image comparison
% imA = MaskASL;
% imB = MaskTemplate;
% TanimotoCoeff = xASL_im_ComputeTanimoto(imA,imB);

if bRegistrationControl
% 1) Initial mean control registrations, if available
    if ~xASL_exist(x.P.Path_mean_control,'file') && xASL_exist(x.P.Path_M0)
        fprintf('No control image present, running M0-T1w registration\n');
        SourcePath = x.P.Path_M0;
        OtherList = xASL_adm_RemoveFromOtherList(BaseOtherList, {x.P.Path_M0});
    elseif ~xASL_exist(x.P.Path_mean_control,'file')
        warning('Skipping control-T1w or M0-T1w registration, couldnt find images');
        SourcePath = NaN;
    else
        fprintf('Running Control-T1w registration\n');
        SourcePath = x.P.Path_mean_control;
        OtherList = xASL_adm_RemoveFromOtherList(BaseOtherList, {x.P.Path_mean_control});
    end

    if min(~isnan(SourcePath)) % if we have a control or M0 image for registration
        xASL_im_BackupAndRestoreAll(BaseOtherList, 1); % First backup all NIfTIs & .mat sidecars of BaseOtherList

        % then register
        if ~x.settings.Quality
            xASL_spm_coreg(x.P.Path_T1, SourcePath, OtherList, x, [9 6]);
        else
            xASL_spm_coreg(x.P.Path_T1, SourcePath, OtherList, x);
        end
        TanimotoPerc(end+1) = xASL_im_GetSpatialOverlapASL(x); % get new overlap score

        if TanimotoPerc(end)>=TanimotoPerc(end-1)
            % if alignment improved or remained more or less the same
            xASL_im_BackupAndRestoreAll(BaseOtherList, 3); % delete backup
        else % if alignment got significantly (>5% Tanimoto) worse
            xASL_im_BackupAndRestoreAll(BaseOtherList, 2); % restore NIfTIs from backup
            TanimotoPerc = TanimotoPerc(1:end-1); % remove last iteration
        end
    end
end



%% ----------------------------------------------------------------------------------------
%% 3.   Registration CBF->pseudoCBF
if bRegistrationCBF

    spatCoVit = xASL_im_GetSpatialCovNativePWI(x);
    if x.modules.asl.bRegistrationContrast==3 || x.modules.asl.bContainsDeltaM
        nIT = 2; % force CBF-pGM
        fprintf('\n%s\n\n','x.modules.asl.bRegistrationContrast==3, forcing CBF-based registration irrespective of sCoV');
    elseif spatCoVit>0.667
        nIT = 0;
        fprintf('%s\n','High spatial CoV, skipping CBF-based registration');
    elseif ~x.settings.Quality
        nIT = 1; % speed up for low quality
    else
        nIT = 2;
	end

	% Create a local instance of the bAffineRegistration (so that any changes in it are not reflected in the x-struct for outside
	bAffineRegistration = x.modules.asl.bAffineRegistration;

    % 2) Repeat CBF registrations, with iteratively better estimate of the
    % vascular/tissue perfusion ratio of the template
    if nIT>0
        % keep this same for all sequences, 3D sequences will simply have a lower spatial CoV because of smoothness
        bSkipThis = false;
		for iT=1:nIT
			if ~bSkipThis

				OtherList = xASL_adm_RemoveFromOtherList(BaseOtherList, {x.P.Path_mean_PWI_Clipped});
				xASL_im_BackupAndRestoreAll(BaseOtherList, 1); % First backup all NIfTIs & .mat sidecars of BaseOtherList

				xASL_im_CreatePseudoCBF(x, spatCoVit(end)); % because this scales the mean_PWI_Clipped, this needs to be run after backing up

				% then register
				xASL_spm_coreg(x.P.Path_PseudoCBF, x.P.Path_mean_PWI_Clipped, OtherList, x);
				% and check for improvement
				TanimotoPerc(end+1) = xASL_im_GetSpatialOverlapASL(x); % get new overlap score

				if x.modules.asl.bRegistrationContrast~=3 % if we don't don't force CBF-pGM registration
					if TanimotoPerc(end)>=TanimotoPerc(end-1)
						% if alignment improved or remained more or less the same
						xASL_im_BackupAndRestoreAll(BaseOtherList, 3); % delete backup
					else
						% if alignment got significantly (>1% Tanimoto) worse
						% we don't force CBF-pGM registration
						xASL_im_BackupAndRestoreAll(BaseOtherList, 2); % restore NIfTIs from backup
						bSkipThis = true; % skip next iteration
						if iT == 1
							bAffineRegistration = 0; % skip affine registration and therefore also DCT - only when it fails to improve on the first, not on the second
						end
						TanimotoPerc = TanimotoPerc(1:end-1); % remove last iteration
					end
				end

				spatCoVit(iT+1) = xASL_im_GetSpatialCovNativePWI(x);
			end
		end

        %% Affine registration
		% Note that this is only done upon request (x.modules.asl.bAffineRegistration, advanced option),
        % hence this doesn't have the automatic backup & restore,
        % as the CBF->pseudoCBF registration has above
		if bAffineRegistration==2 % only do affine for high quality processing & low spatial CoV
			bAffineRegistration = spatCoVit(end)<0.4;
		%   else
			% For bAffineRegistration == 1, do always
			% For bAffineRegsitration == 0, do never
		end

        if bAffineRegistration % perform affine or affine+DCT registration

			if x.modules.asl.bDCTRegistration == 0
				% The affine registration option
				fprintf('%s\n','Performing affine registration');

				xASL_im_BackupAndRestoreAll(BaseOtherList, 1); % First backup all NIfTIs & .mat sidecars of BaseOtherList

				xASL_im_CreatePseudoCBF(x, spatCoVit(end));

				% apply also to mean_PWI_clipped and other files
				xASL_spm_affine(x.P.Path_mean_PWI_Clipped, x.P.Path_PseudoCBF, 5, 5, BaseOtherList);

				% Verify if the affine registration improved the alignment
				TanimotoPerc(end+1) = xASL_im_GetSpatialOverlapASL(x); % get new overlap score
				if TanimotoPerc(end)>=TanimotoPerc(end-1)*0.99
					% if alignment improved or remained more or less the same
					xASL_im_BackupAndRestoreAll(BaseOtherList, 3); % delete backup
				else
					% if alignment got significantly (>1% Tanimoto) worse
					xASL_im_BackupAndRestoreAll(BaseOtherList, 2); % restore NIfTIs from backup
					TanimotoPerc = TanimotoPerc(1:end-1); % remove last iteration
				end
			else
				% The affine+DCT registration option
				fprintf('%s\n','Performing affine+DCT registration');

				% Affine+DCT option does not need a backup because no function is modified, but rather a _sn.mat file
				% is created and can be simply deleted if needed
				if x.modules.asl.bDCTRegistration == 1
					xASL_im_CreatePseudoCBF(x, spatCoVit(end));

					% Use Affine with DCT registration as well
					xASL_spm_affine(x.P.Path_mean_PWI_Clipped, x.P.Path_PseudoCBF, 5,5, [], 1, x.settings.Quality);
				else
					% Use Affine with DCT registration with PVC to prepare the contrast
					% Iterate two times to best use the PVC feature
					for iTDCT = 1:2
						xASL_im_CreatePseudoCBF(x, spatCoVit(end),1);
						xASL_spm_affine(x.P.Path_mean_PWI_Clipped, x.P.Path_PseudoCBF, 5,5, [], 1, x.settings.Quality);
					end
				end

				% Verify if the DCT+affine registration improved the alignment
				TanimotoPerc(end+1) = xASL_im_GetSpatialOverlapASL(x); % get new overlap score
				% No need to delete backup if all went fine.
				if TanimotoPerc(end)<TanimotoPerc(end-1)*0.99
					% if alignment got significantly (>1% Tanimoto) worse
					TanimotoPerc = TanimotoPerc(1:end-1); % remove last iteration

					[FpathSnMat, FfileSnMat] = xASL_fileparts(x.P.Path_mean_PWI_Clipped);
					delete(fullfile(FpathSnMat, [FfileSnMat '_sn.mat']));
				end
			end

            spatCoVit(iT+2) = xASL_im_GetSpatialCovNativePWI(x);
        else
            fprintf('%s\n','Skipping affine registration');
        end
    end
    fprintf('\n%s\n','--------------------------------------------------------------------');

    fprintf('%s\n',[num2str(length(TanimotoPerc)) ' registration iterations:']);
    for iT=1:length(TanimotoPerc)
        fprintf('%s\n',['Iteration ' num2str(iT) ', Tanimoto coefficient = ' num2str(100*TanimotoPerc(iT),3) '%']);
    end
    fprintf('%s\n\n','--------------------------------------------------------------------');
end

% Write the Tanimoto coefficient to the output QC structure
x.Output.ASL.TC_ASL2T1w_Perc = TanimotoPerc(end);

%% ----------------------------------------------------------------------------------------
%% Delete temporary files
if x.settings.DELETETEMP
    File2Del = {x.D.Mean_Native, x.D.Bias_Native, x.D.Vasc_Native, x.D.Mask_Native, x.D.raw_Native, x.P.Path_mean_PWI_Clipped,x.P.Path_mean_PWI_Clipped_DCT,...
        x.P.Path_mean_control, x.P.Path_PseudoCBF, x.D.PathMask, x.D.Path_PseudoTissue, x.D.PathMask, x.P.Path_mean_PWI_Clipped_ORI};
    for iL=1:length(File2Del)
        xASL_delete(File2Del{iL});
    end
end


end



%% ==========================================================================================================
%% ==========================================================================================================
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
        IndexList = 2:length(OtherList);
    elseif IndexIs==length(OtherList)
        IndexList = 1:length(OtherList)-1;
    else
        IndexList = [1:IndexIs-1 IndexIs+1:length(OtherList)];
    end
    OtherList = OtherList(IndexList);
end

end


%% ==========================================================================================================
%% ==========================================================================================================
function [TanimotoCoeff, DiceCoeff] = xASL_im_GetSpatialOverlapASL(x)
%xASL_im_GetSpatialOverlapASL Compute the overlap between two images (using
% TC by default)

if ~isfield(x,'ComputeDiceCoeff')
    x.ComputeDiceCoeff = 0; % the PWI masking doesnt really work
    DiceCoeff = NaN;
end

%% Admin
PathMaskTemplate = fullfile(x.dir.SESSIONDIR, 'Mask_Template.nii');
PathTemplate = fullfile(x.dir.SESSIONDIR, 'Mean_CBF_Template.nii');
[Fpath, Ffile] = xASL_fileparts(x.D.PathMask);
x.D.PathMask2 = fullfile(Fpath, [Ffile '2.nii']);
x.D.PathCBF = fullfile(x.dir.SESSIONDIR, 'LowRes_Mean_CBF_Template.nii');
if ~xASL_exist(x.D.PathMask,'file')
    xASL_Copy(PathMaskTemplate, x.D.PathMask);
end
PWIim = xASL_io_Nifti2Im(x.P.Path_mean_PWI_Clipped);

xASL_spm_reslice(x.P.Path_mean_PWI_Clipped, PathMaskTemplate, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.settings.Quality, x.D.PathMask, 0);
xASL_spm_reslice(x.P.Path_mean_PWI_Clipped, PathTemplate, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.settings.Quality, x.D.PathCBF, 1);

MaskFromTemplate = xASL_io_Nifti2Im(x.D.PathMask)>0.5;
TemplateIm = xASL_io_Nifti2Im(x.D.PathCBF);


%% Compute wholebrain Tanimoto coefficient
TanimotoCoeff = xASL_qc_TanimotoCoeff(PWIim, TemplateIm, MaskFromTemplate, 3, 0.975);
fprintf('%s\n',['Tanimoto Coeff=' num2str(100*TanimotoCoeff,3)]);

if x.ComputeDiceCoeff
    xASL_spm_reslice(x.P.Path_mean_PWI_Clipped, x.D.Mean_Native, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.settings.Quality, x.D.PathMask2 ,0);

    GMIM = xASL_io_Nifti2Im(x.D.PathMask2);
    xASL_delete(x.D.PathMask2);

    %% Compute Dice coefficient
    % Mask the GMIM
    sortInt = sort(GMIM(:));
    ThrInt = sortInt(round(0.7*length(sortInt)));
    GMIM = GMIM>ThrInt;

    % Ensure that the mask is binary
    MaskFromTemplate = MaskFromTemplate & GMIM;

    %% Simple intersection check
    sortInt = sort(PWIim(isfinite(PWIim)));
    ThrInt = sortInt(round(0.75*length(sortInt)));
    maskPWI = PWIim>ThrInt;

    DiceCoeff = xASL_im_ComputeDice(maskPWI,MaskFromTemplate);
    fprintf('%s\n',['Joint brainmask Dice= ' num2str(100*DiceCoeff,3) '%']);
end
xASL_delete(x.D.PathCBF);


end


%% ==========================================================================================================
%% ==========================================================================================================
function [spatCoV] = xASL_im_GetSpatialCovNativePWI(x)
%xASL_im_GetSpatialCovNativePWI Acquires spatial CoV from the native space ASL
%image, using registered mask

JointMasks = xASL_im_GetSpatialOverlapASL(x);

PWIim = xASL_io_Nifti2Im(x.P.Path_mean_PWI_Clipped);
MaskIM = xASL_io_Nifti2Im(x.D.PathMask)>0.5;

if JointMasks<0.5
    warning('Registration off, spatial CoV detection unreliable');
	fprintf('Consider using another registration option\n');
end

%% Determine spatial CoV
spatCoV = xASL_stat_ComputeSpatialCoV(PWIim, MaskIM, 0);

spatCoV = spatCoV/1.5; % correction native space including WM to MNI spatial CoV, excluding WM

if spatCoV<0
    warning('Native space whole-brain spatial CoV was negative! (i.e. <0)');
    fprintf('%s\n', 'Defaulting to spatial CoV of 40%');
    spatCoV = 0.4;
end

fprintf('%s\n', ['Standard space whole-brain spatial CoV estimated as = ' num2str(100*spatCoV,3) '%']);

end



%% ==========================================================================================================
%% ==========================================================================================================
function xASL_im_BackupAndRestoreAll(BaseOtherList, Option)

% Option 1 = backup
% Option 2 = restore from backup
% Option 3 = remove backup

for iNii=1:length(BaseOtherList)
    [Fpath, Ffile, Fext] = xASL_fileparts(BaseOtherList{iNii});
    OriNiiName = BaseOtherList{iNii};
    BackupNiiName = fullfile(Fpath, [Ffile '_Backup' Fext]);
    OriMatName = fullfile(Fpath, [Ffile '.mat']);
    BackupMatName = fullfile(Fpath, [Ffile '_Backup.mat']);

    if Option==1
        % backup
        if xASL_exist(OriNiiName, 'file')
            xASL_Copy(OriNiiName, BackupNiiName, 1, 0);
        end
        if xASL_exist(OriMatName, 'file')
            xASL_Copy(OriMatName, BackupMatName, 1, 0);
        end
    elseif Option==2
        % restore from backup
        if xASL_exist(BackupNiiName, 'file')
            xASL_Move(BackupNiiName, OriNiiName, 1, 0);
        end
        if xASL_exist(BackupMatName, 'file')
            xASL_Move(BackupMatName, OriMatName, 1, 0);
        end
    elseif Option==3
        % remove backup
        xASL_delete(BackupNiiName);
        xASL_delete(BackupMatName);
    end
end


end
