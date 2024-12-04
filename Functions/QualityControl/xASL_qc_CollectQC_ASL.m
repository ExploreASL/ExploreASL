function [x] = xASL_qc_CollectQC_ASL(x, iSubject, iSession)
%xASL_qc_CollectQC_ASL Collect ASL QC parameters
%
% FORMAT: [x] = xASL_qc_CollectQC_ASL(x, iSubject, iSession)
%
% INPUT:
%   x 	     - structure containing fields with all information required to run this submodule (REQUIRED)
%   iSubject - index of current subject (REQUIRED)
%   iSession - index of current session (REQUIRED)
%
% OUTPUT:
%   x        - same as input
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This functions collects QC parameters for the ASL module
%
%              These are stored in x.Output.ASL:
%
%              ID - SubjectName
%              ASL_LR_flip_YesNo - Checks whether any image processing changed the left-right orientation
%                                  by checking whether the determinant differs between nii.mat & nii.mat0
%              SPM realign (too much motion is suspicious)
%               MotionMean_mm    - mean motion
%               MotionExcl_Perc  - percentage of excluded outliers
%               MotionMax_mm     - max motion
%               MotionSD_mm      - SD motion
%
%              ASL quantification (strange average CBF, or strange GM-WM contrast)
%              ASL acquisition parameters (should be fairly consistent over subjects/scans):
%               TE - echo time
%               TR - repetition time
%               RescaleSlope - Philips
%               Scaleslope - Philips
%               Matrix X Y Z - matrix size
%               Matrix Z - number of slices
%               VoxelSize X Y - in plane resolution
%               VoxelSize Z - slice thickness
%               RigidBody2Anat_mm - Net Displacement Vector (RMS) from ASL to T1w image (mm) from registration
%
% With the following parameters:
% 0. Admin
% 1. ASL determinant (left-right flip)
% 2. ASL motion
% 3. Calculate ASL derivatives
% 4. ASL acquisition parameters
% 5. Compute orientation stuff
% 6. RMS, AI, etc of ASL data
% 7. Add data to the QC fields
%
% EXAMPLE: x = xASL_qc_CollectQC_ASL(x, 10, 1);
% __________________________________
% Copyright (c) 2015-2024 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________



    %% 0. Admin
    ASL = struct;
    SubjectID = x.SUBJECTS{iSubject};
    SessionID = x.SESSIONS{iSession};
    ASL_ID = [SubjectID '_' SessionID];


    %% 1. ASL determinant (left-right flip)
    % The determinant of the current matrix and old matrix should be the same,
    % otherwise this is suspicious of a left-right flip.
    PathOrientationResults = fullfile(x.dir.SESSIONDIR,'xASL_qc_PrintOrientation_RigidRegASL.tsv');
    ASL.LR_flip_YesNo = uint8(xASL_im_DetermineFlip(PathOrientationResults));

    if ASL.LR_flip_YesNo>0
        fprintf(['LR flip found for ' SubjectID '_' SessionID]);
    end


    %% 2. ASL motion
    PathMoCo = fullfile(x.D.MotionDir,['motion_correction_NDV_' ASL_ID '.mat']);
    if exist(PathMoCo,'file')
        MoCo = load(PathMoCo);
        ASL.MotionMean_mm = xASL_round(MoCo.mean_NDV{2},4);
        ASL.MotionExcl_Perc = xASL_round(MoCo.PercExcl,3);
        ASL.MotionMax_mm = xASL_round(MoCo.max_NDV{2},4);
        ASL.MotionSD_mm = xASL_round(MoCo.SD_NDV{2},4);
    else
        ASL.MotionMean_mm = NaN;
        ASL.MotionExcl_Perc = NaN;
        ASL.MotionMax_mm = NaN;
        ASL.MotionSD_mm = NaN;
    end


    %% 3. Calculate ASL derivatives
    ASL = xASL_qc_CollectQC_ASL_CalculateDerivatives(x, ASL);
    % Here we calculate a lot of native space ASL derivatives that can be used for QC
    % or other analyses without running the Population module


    %% 4. ASL acquisition parameters
    KnownUnits = {'EchoTime' 'RepetitionTime' 'LabelingDuration' 'Initial_PLD'  'TotalReadoutTime' 'AcquisitionTime' 'SliceReadoutTime'};
    HaveUnits = {'ms'       'ms'             'ms'               'ms'            's'                'hhmmss'          'ms'};

    if isfield(x,'Q')
        QuantFields = fields(x.Q); % all quantification fields
        for iField = 1:length(QuantFields) % iterate over fields
            FieldName = QuantFields{iField};
            IndexIs = find(cellfun(@(x) strcmp(x,FieldName), KnownUnits)); % check if we know the unit
            if ~isempty(IndexIs) % do we know the unit?
                FieldName = [FieldName '_' HaveUnits{IndexIs}]; % then add the unit to the fieldname
            end
            ASL.(FieldName) = x.Q.(QuantFields{iField}); % add the field to ASL struct
        end
    end

    %% 5. Compute orientation stuff
    ASL = xASL_qc_ComputeNiftiOrientation(x.P.Path_ASL4D, ASL);
    
    %% 6. RMS, AI, etc of ASL data
    if strcmp(SessionID(1:3), 'ASL')
        QC_diff_template = xASL_qc_CompareTemplate(x, 'qCBF', iSubject);
    else
        QC_diff_template = xASL_qc_CompareTemplate(x, 'mean_control', iSubject);
    end
    InputFields = fields(QC_diff_template); % add fields to ASL
    for iL=1:length(InputFields)
        if ~isfield(ASL,InputFields{iL})
            ASL.(InputFields{iL}) = QC_diff_template.(InputFields{iL});
        end
    end

    %% 7. Add data to the QC fields
    % Set ASL fields to 4 decimals
    FieldNames = fields(ASL);
    for iN=1:length(FieldNames)
        V = ASL.(FieldNames{iN});
        if isnumeric(V)
            ASL.(FieldNames{iN}) = xASL_round(V, 4);
        end
    end

    % Add data to the QC fields
    Field2Check = fields(ASL);
    nFields = length(Field2Check);
    SumData = 0;
    for iL=1:nFields
        if ~strcmp(Field2Check{iL},'ID') && ~isstruct( ASL.(Field2Check{iL}) )
            if isnumeric( ASL.(Field2Check{iL}) )
                   SumData = SumData+1;
            elseif isnan( ASL.(Field2Check{iL}) )
                   SumData = SumData+1;
            end
        end
	end

	% Check for a session subfield and create when necessary
	if ~isfield(x.Output.ASL, SessionID)
		x.Output.ASL.(SessionID) = struct;
	end

    FieldsFilled = SumData/nFields;
    if FieldsFilled>0.2 % threshold to avoid listing empty values
        x.Output.ASL.(SessionID) = xASL_qc_FillFields(x.Output.ASL.(SessionID), ASL);
    end

end


%% ========================================================================================
%% ========================================================================================


function [OutputFields] = xASL_qc_FillFields(OutputFields, InputFields)
%xASL_qc_FillFields Fill fields        
FieldsI = fields(InputFields);

for iO=1:length(FieldsI)
    CurrentField = InputFields.(FieldsI{iO});
    if isnumeric(CurrentField) && length(CurrentField)>1
        CurrentField = num2str(CurrentField);
    end
    OutputFields.(FieldsI{iO}) = CurrentField;
end    

end


%% ========================================================================================
%% ========================================================================================


function [ASL] = xASL_qc_CollectQC_ASL_CalculateDerivatives(x, ASL)
%xASL_qc_CollectQC_ASL_CalculateDerivatives Calculate ASL parameters
%
% With the following steps:
% 1. Admin
% 2. Calculations over full time-series, excluding vascular signal
% 3. Calculations over full time-series, including vascular signal
% 4. Calculations across time (temporal analyses), including vascular signal
%   I. Admin
%   II. Mask CBF4D but reshape the vector to have the time information in the 2nd dimension
%   III. Tissue masking of the time-series
%   IV. Get spatial parameters per pair
%   V. Calculate temporal values
%   VI. Mean of temporal SD



    %% Get CBF & spatial CoV
    fprintf('%s\n', 'ASL QC: computing native space ASL derivatives...');

    %% 1. Admin
	if xASL_exist(x.P.Path_c1T1,'file') && xASL_exist(x.P.Path_c2T1,'file')
		Path_pGM = x.P.Path_PVgm;
		Path_pWM = x.P.Path_PVwm;
    else
        x = xASL_adm_DefineASLResolution(x);
		warning('T1w-related files missing, computing ASL data using MNI templates!!!');
        xASL_im_PreSmooth(x.P.Path_CBF,fullfile(x.D.TemplateDir,'rc1T1_ASL_res.nii'),...
            x.P.Path_rc1T1,x.S.optimFWHM_Res_mm,[],x.P.Path_mean_PWI_Clipped_sn_mat, 1);
        xASL_im_PreSmooth(x.P.Path_CBF,fullfile(x.D.TemplateDir,'rc2T1_ASL_res.nii'),...
            x.P.Path_rc2T1,x.S.optimFWHM_Res_mm,[],x.P.Path_mean_PWI_Clipped_sn_mat, 1);
	
		xASL_spm_reslice(x.P.Path_CBF, x.P.Path_rc1T1, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.settings.Quality, x.P.Path_rc1T1);
		xASL_spm_reslice(x.P.Path_CBF, x.P.Path_rc2T1, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.settings.Quality, x.P.Path_rc2T1);
		
		Path_pGM = x.P.Path_rc1T1;
		Path_pWM = x.P.Path_rc2T1;
	end
	
    % Read CBF NIfTI
    pGM = xASL_io_Nifti2Im(Path_pGM);
    pWM = xASL_io_Nifti2Im(Path_pWM);
    imCBF = xASL_io_Nifti2Im(x.P.Path_CBF);
	imCBF4D = xASL_io_Nifti2Im(x.P.Path_CBF4D);
        
    if xASL_stat_SumNan(pGM(:))==0
        warning(['Empty image, invalid ' Path_pGM]);
    end
    if xASL_stat_SumNan(pWM(:))==0
        warning(['Empty image, invalid ' Path_pWM]);
    end
    if xASL_stat_SumNan(imCBF(:))==0
        warning(['Empty image, invalid ' x.P.Path_CBF]);
	end  
	if xASL_stat_SumNan(imCBF4D(:))==0
        warning(['Empty image, invalid ' x.P.Path_CBF4D]);
    end    
    

    %% 2. Calculations over full time-series, excluding vascular signal
    % NB: x.P.Path_MaskVascular needs to be present!

    if xASL_exist(x.P.Path_MaskVascular, 'file')
        imMaskWB = (pGM+pWM)>0.5;
        imMaskWB = logical(imMaskWB.*(xASL_io_Nifti2Im(x.P.Path_MaskVascular)>0));
        CBFmasked = imCBF(imMaskWB);
        GMmasked = pGM(imMaskWB);
        WMmasked = pWM(imMaskWB);        
    
        % CBF
        ASL.CBF_GM_Median_mL100gmin = xASL_stat_ComputeMean(CBFmasked, GMmasked>0.7,[], 0, 0);
        % PM: this name should be changed later

        [ASL.CBF_GM_PVC2_mL100gmin, ASL.CBF_WM_PVC2_mL100gmin] = xASL_stat_ComputeMean(CBFmasked, (GMmasked+WMmasked)>0.5,[],2, 1, GMmasked, WMmasked);
        ASL.CBF_GM_WM_Ratio = ASL.CBF_GM_PVC2_mL100gmin/ASL.CBF_WM_PVC2_mL100gmin;
    else
        warning(['Missing: ' x.P.Path_MaskVascular]);
        fprintf('%s\n', 'Need vascular mask to calculate native space CBF values!');
    end

    
    %% 3. Calculations over full time-series, including vascular signal
    
    % Spatial CoV
    CBFmasked = imCBF(imMaskWB);
    ASL.SpatialCoV_GM_Perc = 100*xASL_stat_ComputeSpatialCoV(CBFmasked, [], [], 0, 1);


    %% 4. Calculations across time, including vascular signal
    %% I. Admin
    
    nPairs = size(imCBF4D, 4);
    
    imMaskWB = (pGM+pWM)>0.5;
    WBmasked = logical(ones([sum(imMaskWB(:)), 1]));
    GMmasked = pGM(imMaskWB) > 0.7; % same as used above
    WMmasked = pWM(imMaskWB) > 0.7; % fits approx. with pGM>0.5

    WBmasked4D = repmat(WBmasked, [1 nPairs]);
    GMmasked4D = repmat(GMmasked, [1 nPairs]);
    WMmasked4D = repmat(WMmasked, [1 nPairs]);

    %% II. Mask CBF4D but reshape the vector to have the time information in the 2nd dimension
	% Get all voxels within the mask
    CBF4Dmasked = imCBF4D(repmat(imMaskWB, [1 1 1 nPairs]));
    % Reshape the vector
    CBF4Dmasked = reshape(CBF4Dmasked, [], nPairs);

	%% III. Tissue masking of the time-series
	CBF4DmaskedWB = reshape(CBF4Dmasked(WBmasked4D), [], nPairs);
	CBF4DmaskedGM = reshape(CBF4Dmasked(GMmasked4D), [], nPairs);
	CBF4DmaskedWM = reshape(CBF4Dmasked(WMmasked4D), [], nPairs);

    %% IV. Get spatial parameters per pair
    
    % PM: For now restricted to whole-brain WB only, for simplicity
    % As we need to include both the CoW and distal areas

	CBF.mean = xASL_stat_MeanNan(CBF4DmaskedWB, 1);
	CBF.meanGM = xASL_stat_MeanNan(CBF4DmaskedGM, 1);
	CBF.median = xASL_stat_MedianNan(CBF4DmaskedWB, 1);
	CBF.SD = xASL_stat_StdNan(CBF4DmaskedWB, [], 1);
	CBF.MAD = xASL_stat_MadNan(CBF4DmaskedWB, 1);
	CBF.sCoV = CBF.SD./CBF.mean;
	for iRepetition=1:nPairs
        CBF.diffCoV(iRepetition) = xASL_stat_ComputeDifferCoV(imCBF4D(:, :, :, iRepetition), imMaskWB);
	end

    %% V. Calculate temporal values

    % PM: keeping it simple here for now, we can add more parameters later
    % Naming: 
    % ASL.<parameter calculated first>_<over which ROI>_temporal<parameter calculated second>

    ASL.SpatialCoV_WB_temporalMean = 100.*xASL_stat_MeanNan(CBF.sCoV);
    ASL.SpatialCoV_WB_temporalSD = 100.*xASL_stat_StdNan(CBF.sCoV);

    ASL.SpatialSD_WB_temporalMean = xASL_stat_MeanNan(CBF.SD);
    ASL.SpatialSD_WB_temporalSD = xASL_stat_StdNan(CBF.SD);

    ASL.DiffCoV_WB_temporalMean = 100.*xASL_stat_MeanNan(CBF.diffCoV);
    ASL.DiffCoV_WB_temporalSD = 100.*xASL_stat_StdNan(CBF.diffCoV);
	
	%% VI. Mean of temporal SD
	ASL.tSD_WB_Mean = xASL_stat_MeanNan(xASL_stat_StdNan(CBF4DmaskedWB, [], 2), 1);
    ASL.tSD_GM_Mean = xASL_stat_MeanNan(xASL_stat_StdNan(CBF4DmaskedGM, [], 2), 1);
	ASL.tSD_WM_Mean = xASL_stat_MeanNan(xASL_stat_StdNan(CBF4DmaskedWM, [], 2), 1);


end