function [x] = xASL_qc_CollectQC_func(x, iSubject)
%xASL_qc_CollectQC_func Collect func QC parameters
%
% FORMAT: [x] = xASL_qc_CollectQC_func(x, iSubject)
%
% INPUT:
%   x 	     - structure containing fields with all information required to run this submodule (REQUIRED)
%   iSubject - index of current subject (REQUIRED)
%
% OUTPUT:
%   x        - same as input
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This functions collects QC parameters for the func module
%
%              These are stored in x.Output.func:
%
%              ID - SubjectName
%              func_LR_flip_YesNo - Checks whether any image processing changed the left-right orientation
%                                  by checking whether the determinant differs between nii.mat & nii.mat0
%              SPM realign (too much motion is suspicious)
%               MotionMean_mm    - mean motion
%               MotionExcl_Perc  - percentage of excluded outliers
%               MotionMax_mm     - max motion
%               MotionSD_mm      - SD motion
%              func quantification (strange average CBF, or strange GM-WM contrast)
%               CBF_GM_Median_mL100gmin - median GM CBF
%               CBF_WM_Median_mL100gmin - median WM CBF
%               SpatialCoV_GM_Perc      - GM spatial CoV
%               SpatialCoV_WM_Perc      - WM spatial CoV
%               CBF_GM_WM_Ratio         - GM-WM CBF ratio
%
%              func acquisition parameters (should be fairly consistent over subjects/scans):
%               TE - echo time
%               TR - repetition time
%               RescaleSlope - Philips
%               Scaleslope - Philips
%               Matrix X Y Z - matrix size
%               Matrix Z - number of slices
%               VoxelSize X Y - in plane resolution
%               VoxelSize Z - slice thickness
%               RigidBody2Anat_mm - Net Displacement Vector (RMS) from func to T1w image (mm) from registration
%
% EXAMPLE: x = xASL_qc_CollectQC_func(x, 10);
% __________________________________
% Copyright (C) 2015-2019 Explorefunc


    %% Admin
    SubjectID = x.SUBJECTS{iSubject};
    SessionID = x.SESSIONS{1};
    func.ID = [SubjectID '_' SessionID];


    %% -----------------------------------------------------------------------------------------------
    %% func determinant
    % The determinant of the current matrix and old matrix should be the same,
    % otherwise this is suspicious of a left-right flip.

    PathOrientationResults = fullfile(x.SESSIONDIR,'xASL_qc_PrintOrientation_RigidRegfunc.tsv');
    func = xASL_im_DetermineFlip(x, iSubject, PathOrientationResults, func);

    %% func motion
    PathMoCo = fullfile(x.D.MotionDir,['motion_correction_NDV_' func.ID '.mat']);
    if  exist(PathMoCo,'file')
        MoCo                    = load(PathMoCo);
        func.MotionMean_mm       = xASL_round(MoCo.mean_NDV{2},4);
        func.MotionExcl_Perc     = xASL_round(MoCo.PercExcl,3);
        func.MotionMax_mm = xASL_round(MoCo.max_NDV{2},4);
        func.MotionSD_mm = xASL_round(MoCo.SD_NDV{2},4);        
    else
        func.MotionMean_mm       = NaN;
        func.MotionExcl_Perc     = NaN;
        func.MotionMax_mm = NaN;
        func.MotionSD_mm = NaN;        
    end


    %% -----------------------------------------------------------------------------------------------
    %% func acquisition
    
    KnownUnits = {'EchoTime' 'RepetitionTime' 'TotalReadoutTime' 'AcquisitionTime'};
    HaveUnits = {'ms'       'ms'             's'                'hhmmss'};

    for iField=1:length(KnownUnits)
        if isfield(x, KnownUnits{iField}) && ~isfield(x.Q, KnownUnits{iField})
            x.Q.(KnownUnits{iField}) = x.(KnownUnits{iField});
        end
    end
    
    if isfield(x,'Q')
        QuantFields = fields(x.Q); % all quantification fields
        for iField = 1:length(QuantFields) % iterate over fields
            FieldName = QuantFields{iField};
            IndexIs = find(cellfun(@(x) strcmp(x,FieldName), KnownUnits)); % check if we know the unit
            if ~isempty(IndexIs) % do we know the unit?
                FieldName = [FieldName '_' HaveUnits{IndexIs}]; % then add the unit to the fieldname
            end
            func.(FieldName) = x.Q.(QuantFields{iField}); % add the field to ASL struct
        end
    end    

    % Orientation check
    func = xASL_qc_ComputeNiftiOrientation(x, x.P.Path_func_bold, func);
    
    %% Set func fields to 4 decimals
    FieldNames = fields(func);
    for iN=1:length(FieldNames)
        V = func.(FieldNames{iN});
        if isnumeric(V)
            func.(FieldNames{iN}) = xASL_round(V, 4);
        end
    end

    % First check if we have some data to add
    Field2Check     = fields(func);
    nFields         = length(Field2Check);
    SumData         = 0;
    for iL=1:nFields
        if ~strcmp(Field2Check{iL},'ID') && ~isstruct( func.(Field2Check{iL}) )
            if      isnumeric( func.(Field2Check{iL}) )
                    SumData     = SumData+1;
            elseif  isnan( func.(Field2Check{iL}) )
                    SumData     = SumData+1;
            end
        end
    end
    FieldsFilled    = SumData/nFields;
    if  FieldsFilled>0.2 % threshold to avoid listing empty values
        x.Output.('func') = xASL_qc_FillFields(x.Output.('func'), func);
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