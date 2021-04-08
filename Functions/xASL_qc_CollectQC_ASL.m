function [x] = xASL_qc_CollectQC_ASL(x, iSubject)
%xASL_qc_CollectQC_ASL Collect ASL QC parameters
%
% FORMAT: [x] = xASL_qc_CollectQC_ASL(x, iSubject)
%
% INPUT:
%   x 	     - structure containing fields with all information required to run this submodule (REQUIRED)
%   iSubject - index of current subject (REQUIRED)
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
% EXAMPLE: x = xASL_qc_CollectQC_ASL(x, 10);
% __________________________________
% Copyright (C) 2015-2019 ExploreASL


    %% Admin
    ASL = struct;

    SubjectID = x.SUBJECTS{iSubject};
    SessionID = x.SESSIONS{1};
    ASL_ID = [SubjectID '_' SessionID];

    %% -----------------------------------------------------------------------------------------------
    %% ASL determinant
    % The determinant of the current matrix and old matrix should be the same,
    % otherwise this is suspicious of a left-right flip.
    PathOrientationResults = fullfile(x.SESSIONDIR,'xASL_qc_PrintOrientation_RigidRegASL.tsv');
    ASL = xASL_im_DetermineFlip(x, iSubject, PathOrientationResults, ASL);

    %% ASL motion
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


    %% -----------------------------------------------------------------------------------------------
    %% ASL CBF values
    %% Get CBF & spatial CoV
    
	if xASL_exist(x.P.Path_c1T1,'file') && xASL_exist(x.P.Path_c2T1,'file')
		Path_pGM = x.P.Path_PVgm;
		Path_pWM = x.P.Path_PVwm;
    else
        x = xASL_adm_DefineASLResolution(x);
		warning('T1w-related files missing, computing ASL data using MNI templates!!!');
        xASL_im_PreSmooth(x.P.Path_CBF,fullfile(x.D.TemplateDir,'rc1T1_ASL_res.nii'),x.P.Path_rc1T1,x.S.optimFWHM_Res_mm,[],x.P.Path_mean_PWI_Clipped_sn_mat, 1);
        xASL_im_PreSmooth(x.P.Path_CBF,fullfile(x.D.TemplateDir,'rc2T1_ASL_res.nii'),x.P.Path_rc2T1,x.S.optimFWHM_Res_mm,[],x.P.Path_mean_PWI_Clipped_sn_mat, 1);
	
		xASL_spm_reslice(x.P.Path_CBF, x.P.Path_rc1T1, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.Quality, x.P.Path_rc1T1);
		xASL_spm_reslice(x.P.Path_CBF, x.P.Path_rc2T1, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.Quality, x.P.Path_rc2T1);
		
		Path_pGM = x.P.Path_rc1T1;
		Path_pWM = x.P.Path_rc2T1;
	end
	
    pGM = xASL_io_Nifti2Im(Path_pGM);
    pWM = xASL_io_Nifti2Im(Path_pWM);
    imCBF = xASL_io_Nifti2Im(x.P.Path_CBF);

    if xASL_stat_SumNan(pGM(:))==0
        warning(['Empty image, invalid ' Path_pGM]);
    end
    if xASL_stat_SumNan(pWM(:))==0
        warning(['Empty image, invalid ' Path_pWM]);
    end
    if xASL_stat_SumNan(imCBF(:))==0
        warning(['Empty image, invalid ' x.P.Path_CBF]);
    end    
    
    imMask = (pGM+pWM)>0.5;
    CBFmasked = imCBF(imMask);
    GMmasked = pGM(imMask);
    WMmasked = pWM(imMask);

    % Including vascular signal
    ASL.SpatialCoV_GM_Perc = 100*xASL_stat_ComputeSpatialCoV(CBFmasked,[],[],0);

    if xASL_exist(x.P.Path_MaskVascular, 'file')
        imMask = logical(imMask.*xASL_io_Nifti2Im(x.P.Path_MaskVascular));
    end
    CBFmasked = imCBF(imMask);
    GMmasked = pGM(imMask);
    WMmasked = pWM(imMask);

    % Excluding vascular signal
    ASL.CBF_GM_Median_mL100gmin = xASL_stat_ComputeMean(CBFmasked, GMmasked>0.5,[], 0, 0);
    [ASL.CBF_GM_PVC2_mL100gmin, ASL.CBF_WM_PVC2_mL100gmin] = xASL_stat_ComputeMean(CBFmasked, (GMmasked+WMmasked)>0.5,[],2, 1, GMmasked, WMmasked);
    ASL.CBF_GM_WM_Ratio = ASL.CBF_GM_PVC2_mL100gmin/ASL.CBF_WM_PVC2_mL100gmin;


    %% -----------------------------------------------------------------------------------------------
    %% ASL acquisition
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

    % compute orientation stuff
    ASL = xASL_qc_ComputeNiftiOrientation(x.P.Path_ASL4D, ASL);
    
    %% RMS, AI, etc of ASL data
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

    %% Set ASL fields to 4 decimals
    FieldNames = fields(ASL);
    for iN=1:length(FieldNames)
        V = ASL.(FieldNames{iN});
        if isnumeric(V)
            ASL.(FieldNames{iN}) = xASL_round(V, 4);
        end
    end

    %% Add data to the QC fields
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
    FieldsFilled = SumData/nFields;
    if FieldsFilled>0.2 % threshold to avoid listing empty values
        x.Output.('ASL') = xASL_qc_FillFields(x.Output.('ASL'), ASL);
    end

end

function [OutputFields] = xASL_qc_FillFields(OutputFields, InputFields)
%xASL_qc_FillFields
        
FieldsI = fields(InputFields);

for iO=1:length(FieldsI)
    CurrentField = InputFields.(FieldsI{iO});
    if isnumeric(CurrentField) && length(CurrentField)>1
        CurrentField = num2str(CurrentField);
    end
    OutputFields.(FieldsI{iO}) = CurrentField;
end    

end