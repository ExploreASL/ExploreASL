function xASL_wrp_RegisterASL(x)
%xASL_wrp_RegisterASL Submodule of ExploreASL ASL Module, that registers
%ASL to T1w (or potentially other structural images)
%
% FORMAT: xASL_wrp_RegisterASL(x)
%
% INPUT:
%   x  - structure containing fields with all information required to run this submodule (REQUIRED)
%
% OUTPUT: n/a (registration changes the NIfTI orientation header only
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
%                           - 0 = Control->T1w

%                           - 1 = CBF->pseudoCBF from template/pGM+pWM
%                                 (skip if sCoV>0.667)
%                           - 2 = automatic (mix of both)
%                           - 3 = option 2 & force CBF->pseudoCBF irrespective of sCoV or Tanimoto coefficient
%     - G) Dummy src NIfTIs are created:
%          mean_control.nii to register with T1w
%          mean_PWI_Clipped.nii to register with pseudoCBF
%     - H) Create reference images, downsampled pseudoTissue
%
% 1)    Registration Center of Mass
% 2)    Registration ASL -> anat (Control->T1w)
%       (in case of a 3D sequence, this step is only applied if it improves the Tanimoto coefficient by more than 1%)
% 3)    Registration CBF->pseudoCBF
%       (in case of a 2D sequence, this step is only applied if it improves
%       the Tanimoto coefficient by more than 1%). Also, this step is only
%       applied if the spatial CoV<0.67. Note that this is usually the case
%       for 3D scans because of their lower effective spatial resolution.
%
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

if ~isfield(x,'bRegistrationContrast')
    x.bRegistrationContrast = 2; % register M0-T1w first, then CBF-pGM if sCoV<0.667
end

%% B) Manage OtherList
% Define OtherList for registration
% Here we mention all possible files that need to be in registration
% All functions below will remove those that are unexisting, or used in the
% registration estimation.
BaseOtherList = {x.P.Path_despiked_ASL4D, x.P.Path_mean_control, x.P.Path_M0, x.P.Path_PWI, x.P.Path_mean_PWI_Clipped, x.P.Path_ASL4D_RevPE,...
    x.P.Path_ASL4D_ORI, fullfile(x.SESSIONDIR, 'B0.nii'), fullfile(x.SESSIONDIR, 'Unwarped.nii'), fullfile(x.SESSIONDIR, 'Field.nii'), fullfile(x.SESSIONDIR, 'TopUp_fieldcoef.nii')};

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
 
    

%% E) Allow registration without structural data
StructuralDerivativesExist = xASL_exist(x.P.Path_y_T1, 'file') && xASL_exist(x.P.Path_c1T1, 'file') && xASL_exist(x.P.Path_c2T1, 'file');
if xASL_exist(x.P.Path_T1, 'file') && ~StructuralDerivativesExist
    error('Please run structural module first');
elseif ~xASL_exist(x.P.Path_T1, 'file') && ~StructuralDerivativesExist
    warning('Missing structural scans, using ASL registration only instead!');
    IDmatrixPath = fullfile(x.D.MapsSPMmodifiedDir, 'Identity_Deformation_y_T1.nii');
    % Copy dummy transformation field
    xASL_Copy(IDmatrixPath, x.P.Path_y_T1, true);
    % Create dummy native space structural derivatives
    % In standard space
	xASL_Copy(fullfile(x.D.MapsSPMmodifiedDir, 'rc1T1.nii'), x.P.Pop_Path_rc1T1);
	xASL_Copy(fullfile(x.D.MapsSPMmodifiedDir, 'rc2T1.nii'), x.P.Pop_Path_rc2T1);
	xASL_Copy(fullfile(x.D.MapsSPMmodifiedDir, 'rT1.nii'), x.P.Pop_Path_rT1);

    % In native space
    xASL_spm_deformations(x, {x.P.Pop_Path_rc1T1, x.P.Pop_Path_rc2T1, x.P.Pop_Path_rT1}, {x.P.Path_c1T1, x.P.Path_c2T1, x.P.Path_T1});
	
	% Dummy files
	catVolFile = fullfile(x.D.TissueVolumeDir,['cat_' x.P.STRUCT '_' x.P.SubjectID '.mat']);
	MatFile   = fullfile(x.SUBJECTDIR, [x.P.STRUCT '_seg8.mat']);
	dummyVar = [];
	save(catVolFile,'dummyVar');
	save(MatFile,'dummyVar');
	
	SaveFile = fullfile(x.D.TissueVolumeDir,['TissueVolume_' x.P.SubjectID '.csv']);
    FileID = fopen(SaveFile,'wt');
	fprintf(FileID,'%s', '0, 0, 0');
	fclose(FileID);
	xASL_adm_csv2tsv(SaveFile, true);
	
	% To lock in the structural part	
    % Save the ASL lock and unlock
	jj = strfind(x.LOCKDIR,'xASL_module_ASL');
	jj = jj(1);
	oldRoot = x.mutex.Root;
	oldID = x.mutex.ID;
	newRoot = fullfile(x.LOCKDIR(1:(jj-1)),'xASL_module_Structural',x.LOCKDIR((jj+16):end));
	x.mutex.Unlock();
	
	% Look the structural part
	x.mutex.Root = newRoot;
	x.mutex.Lock('xASL_module_Structural');

	% Add the correct lock-files
	x.mutex.AddState('010_LinearReg_T1w2MNI');
	x.mutex.AddState('060_Segment_T1w');
	x.mutex.AddState('080_Resample2StandardSpace');
	x.mutex.AddState('090_GetVolumetrics');
	x.mutex.AddState('100_VisualQC_Structural');
	x.mutex.AddState('110_DoWADQCDC');
	
	% Unlock the structural and lock again the ASL part
	x.mutex.Unlock();
	x.mutex.Root = oldRoot;
	x.mutex.Lock(oldID);

end
    
%% E) % Smooth T1 deformation field into ASL resolution
% If no T1 flow field exists, create an identity flowfield
% So we can still process ASL images without the
% structural module
if ~xASL_exist(x.P.Path_y_ASL,'file') || strcmp(x.P.SessionID,'ASL_1') || x.nSessions==1
    xASL_wrp_CreateASLDeformationField(x, true);
end


%% F) Manage registration contrasts that we will use
if x.bRegistrationContrast==0
    bRegistrationControl = true;
    bRegistrationCBF = false;
elseif x.bRegistrationContrast==1
    bRegistrationControl = false;
    bRegistrationCBF = true;
else
    bRegistrationControl = true;
    bRegistrationCBF = true;
end

%% G) Here we create a temporary dummy ASL image of which the image contrast is curated,
% for registration only
xASL_io_PairwiseSubtraction(x.P.Path_despiked_ASL4D, x.P.Path_mean_PWI_Clipped, 0, 0); % create PWI & mean_control
if bRegistrationCBF
    % Clip & compress the image, deal with contrast used for registration
    tIM = xASL_im_ClipExtremes(x.P.Path_mean_PWI_Clipped, 0.95, 0.6); % careful, this cannot be rerun, once
    tIM(tIM==min(tIM(:))) = 0; % minimal intensities are set to 0
    xASL_io_SaveNifti(x.P.Path_mean_PWI_Clipped, x.P.Path_mean_PWI_Clipped, tIM, [], 0);
end

%% H) Here we create the reference images, the downsampled pseudoTissue
% (pGM+pWM) as well as the native space copies of templates for CBF,
% ATT biasfield and vascular peaks
xASL_im_CreatePseudoCBF(x, 0);

TanimotoPerc{1} = xASL_im_GetSpatialOverlapASL(x);


%% ----------------------------------------------------------------------------------------
%% 1)   Registration CenterOfMass
% We now assume the structural module hasn't run, and we simply want to run the ASL module to quickly check how the images look like
% So we only run the automatic Center of Mass ACPC alignment
if x.bAutoACPC
    OtherList = xASL_adm_RemoveFromOtherList(BaseOtherList, {x.P.Path_despiked_ASL4D}); % x.P.Path_despiked_ASL4D is padded to the end of the list
    
    xASL_im_BackupAndRestoreAll(BaseOtherList, 1); % First backup all NIfTIs & .mat sidecars of BaseOtherList
    xASL_im_CenterOfMass(x.P.Path_despiked_ASL4D, OtherList, 0); % Then register
    TanimotoPerc{end+1} = xASL_im_GetSpatialOverlapASL(x); % get new overlap score
    if TanimotoPerc{end}>=TanimotoPerc{end-1} % if alignment improved or remained same
        xASL_im_BackupAndRestoreAll(BaseOtherList, 3); % delete backup
    else % if alignment got worse
        xASL_im_BackupAndRestoreAll(BaseOtherList, 2); % restore NIfTIs from backup
    end
end


%% ----------------------------------------------------------------------------------------
%% 2)    Registration Control->T1w
% Here we first create a mask
% First check the initial alignment, otherwise first register with template
% xASL_spm_reslice(Mask_Native, x.P.Path_mean_PWI_Clipped, x.P.Path_mean_PWI_Clipped_sn_mat, 0, x.Quality, x.P.Path_rmean_PWI_Clipped, 1);
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
    if ~xASL_exist(x.P.Path_mean_control,'file')
        warning('Trying control-T1w registration but Mean_control.nii didnt exist');
    else
        OtherList = xASL_adm_RemoveFromOtherList(BaseOtherList, {x.P.Path_mean_control});

        if isempty(regexp(x.readout_dim,'2D')) % for 2D (EPI) we don't care (assuming this works better than center of mass)
            xASL_im_BackupAndRestoreAll(BaseOtherList, 1); % First backup all NIfTIs & .mat sidecars of BaseOtherList
        end
        
        % then register
        if ~x.Quality
            xASL_spm_coreg(x.P.Path_T1, x.P.Path_mean_control, OtherList, x,[9 6]);
        else
            xASL_spm_coreg(x.P.Path_T1, x.P.Path_mean_control, OtherList, x);
        end
        TanimotoPerc{end+1} = xASL_im_GetSpatialOverlapASL(x); % get new overlap score
        
        if isempty(regexp(x.readout_dim,'2D')) % for 2D (EPI) we don't care (assuming this works better than center of mass)
            if TanimotoPerc{end}>=TanimotoPerc{end-1}
                % if alignment improved or remained more or less the same
                xASL_im_BackupAndRestoreAll(BaseOtherList, 3); % delete backup
                FinalTanimoto = TanimotoPerc{end};
            else % if alignment got significantly (>5% Tanimoto) worse
                xASL_im_BackupAndRestoreAll(BaseOtherList, 2); % restore NIfTIs from backup
                FinalTanimoto = TanimotoPerc{end-1};
            end
        end
    end
end



%% ----------------------------------------------------------------------------------------
%% 3)   Registration CBF->pseudoCBF
if bRegistrationCBF

    spatCoVit = xASL_im_GetSpatialCovNativePWI(x);
    if x.bRegistrationContrast==3
        nIT = 2; % force CBF-pGM
        fprintf('\n%s\n\n','x.bRegistrationContrast==3, forcing CBF-based registration irrespective of sCoV');
    elseif spatCoVit>0.667
        nIT = 0;
        fprintf('%s\n','High spatial CoV, skipping CBF-based registration');
    elseif ~x.Quality
        nIT = 1; % speed up for low quality
    else
        nIT = 2;
    end


    % 2) Repeat CBF registrations, with iteratively better estimate of the
    % vascular/tissue perfusion ratio of the template
    if nIT>0
        % keep this same for all sequences, 3D sequences will simply have a lower spatial CoV because of smoothness
        bSkipThis = false;
        for iT=1:nIT
            if ~bSkipThis
                xASL_im_CreatePseudoCBF(x, spatCoVit(1));
                OtherList = xASL_adm_RemoveFromOtherList(BaseOtherList, {x.P.Path_mean_PWI_Clipped});

                if isempty(regexp(x.Sequence,'spiral')) % if we don't have a 3D spiral sequence
                    xASL_im_BackupAndRestoreAll(BaseOtherList, 1); % First backup all NIfTIs & .mat sidecars of BaseOtherList
                end
                
                % then register
                xASL_spm_coreg(x.P.Path_PseudoCBF, x.P.Path_mean_PWI_Clipped, OtherList, x);
                % and check for improvement
                TanimotoPerc{end+1} = xASL_im_GetSpatialOverlapASL(x); % get new overlap score
                
                if isempty(regexp(x.Sequence,'spiral')) && x.bRegistrationContrast~=3
                    % if we don't have a 3D spiral sequence & don't force CBF-pGM registration
                    if TanimotoPerc{end}>=TanimotoPerc{end-1}
                        % if alignment improved or remained more or less the same
                        xASL_im_BackupAndRestoreAll(BaseOtherList, 3); % delete backup
                        FinalTanimoto = TanimotoPerc{end};
                    else
                        % if alignment got significantly (>1% Tanimoto) worse
                        % we don't force CBF-pGM registration
                        xASL_im_BackupAndRestoreAll(BaseOtherList, 2); % restore NIfTIs from backup
                        bSkipThis = true; % skip next iteration
                        x.bAffineRegistration = 0; % skip affine registration
                        FinalTanimoto = TanimotoPerc{end-1};
                    end
                end

                spatCoVit(iT+1) = xASL_im_GetSpatialCovNativePWI(x);
            end
		end

        %% Affine registration
        if isfield(x,'bAffineRegistration') && ~isempty(x.bAffineRegistration)
            if x.bAffineRegistration==2 % only do affine for high quality processing & low spatial CoV
                bAffineRegistration = x.Quality && spatCoVit(end)<0.4;
            else
                bAffineRegistration = x.bAffineRegistration;
            end
        else
            bAffineRegistration = false; % default
            % default is no affine registration, as the SPM affine COST
            % function is tricky without proper rescaling & testing this
        end

        if bAffineRegistration % perform affine registration
            % NB: here the Tanimoto coefficient check may not work, when e.g.
            % one image is enlarged, it has more overlap
            fprintf('%s\n','Performing affine registration');
            xASL_im_CreatePseudoCBF(x, spatCoVit(end));
            % apply also to mean_PWI_clipped and other files
            xASL_spm_affine(x.P.Path_mean_PWI_Clipped, x.P.Path_PseudoCBF, 5, 5, BaseOtherList);
            spatCoVit(iT+2) = xASL_im_GetSpatialCovNativePWI(x);
        else
            fprintf('%s\n','Skipping affine registration');
        end
    end
    %%%% Another step with SPM_NORMALISE deformations
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
    File2Del = {x.Mean_Native x.Bias_Native x.Vasc_Native x.Mask_Native x.raw_Native x.P.Path_mean_PWI_Clipped x.P.Path_mean_control x.P.Path_PseudoCBF x.PathMask x.Path_PseudoTissue};
    for iL=1:length(File2Del)
        xASL_delete(File2Del{iL});
    end
    xASL_delete(x.PathMask);
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


function [TanimotoCoeff, DiceCoeff] = xASL_im_GetSpatialOverlapASL(x)
 
 
if ~isfield(x,'ComputeDiceCoeff')
    x.ComputeDiceCoeff = 0; % the PWI masking doesnt really work
    DiceCoeff = NaN;
end
 
%% Admin
PathMaskTemplate = fullfile(x.SESSIONDIR, 'Mask_Template.nii');
PathTemplate = fullfile(x.SESSIONDIR, 'Mean_CBF_Template.nii');
[Fpath, Ffile] = xASL_fileparts(x.PathMask);
x.PathMask2 = fullfile(Fpath, [Ffile '2.nii']);
x.PathCBF = fullfile(x.SESSIONDIR, 'LowRes_Mean_CBF_Template.nii');
if ~xASL_exist(x.PathMask,'file')
    xASL_Copy(PathMaskTemplate, x.PathMask);
end
PWIim = xASL_io_Nifti2Im(x.P.Path_mean_PWI_Clipped);
 
xASL_spm_reslice(x.P.Path_mean_PWI_Clipped, PathMaskTemplate, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.Quality, x.PathMask, 0);
xASL_spm_reslice(x.P.Path_mean_PWI_Clipped, PathTemplate, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.Quality, x.PathCBF, 1);
 
MaskFromTemplate = xASL_io_Nifti2Im(x.PathMask)>0.5;
TemplateIm = xASL_io_Nifti2Im(x.PathCBF);
 
 
%% Compute wholebrain Tanimoto coefficient
TanimotoCoeff = xASL_qc_TanimotoCoeff(PWIim, TemplateIm, MaskFromTemplate, 3, 0.975);
fprintf('%s\n',['Tanimoto Coeff=' num2str(100*TanimotoCoeff,3)]);
 
if x.ComputeDiceCoeff
    xASL_spm_reslice(x.P.Path_mean_PWI_Clipped, x.Mean_Native, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.Quality, x.PathMask2 ,0);
    
    GMIM = xASL_io_Nifti2Im(x.PathMask2);
    xASL_delete(x.PathMask2);
    
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



end


function [spatCoV] = xASL_im_GetSpatialCovNativePWI(x)
%xASL_im_GetSpatialCovNativePWI Acquires spatial CoV from the native space ASL
%image, using registered mask

JointMasks = xASL_im_GetSpatialOverlapASL(x);

PWIim = xASL_io_Nifti2Im(x.P.Path_mean_PWI_Clipped);
MaskIM = xASL_io_Nifti2Im(x.PathMask)>0.5;

if JointMasks<0.25
    warning('Registration off, spatial CoV detection unreliable');
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