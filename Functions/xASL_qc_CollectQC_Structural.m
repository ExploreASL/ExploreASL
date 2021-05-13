function [x] = xASL_qc_CollectQC_Structural(x, iSubject)
%xASL_qc_CollectQC_Structural Collect structural/anatomical QC parameters
%
% FORMAT: [x] = xASL_qc_CollectQC_Structural(x, iSubject)
%
% INPUT:
%   x 	     - structure containing fields with all information required to run this submodule (REQUIRED)
%   iSubject - index of current subject (REQUIRED)
%
% OUTPUT:
%   x        - same as input
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This functions collects QC parameters for the structural module
%              These are stored in x.Output.Structural:
%              ID - SubjectName
%              T1w_LR_flip_YesNo - Checks whether any image processing changed the left-right orientation
%                                  by checking whether the determinant differs between nii.mat & nii.mat0
%              LST output:
%               WMH_vol_mL        - WMH volume
%               WMH_n             - WMH number of lesions
%              CAT12 output:
%               T1w_IQR_Perc      - CAT12 IQR quality parameter for T1w
%               volumetric: GM_vol_mL, WM_vol_mL, CSF_vol_mL, ICV_vol_mL, GM_ICV_Ratio
% 
% EXAMPLE: x = xASL_qc_CollectQC_Structural(x, 10);
% __________________________________
% Copyright (C) 2015-2019 ExploreASL


    %% -----------------------------------------------------------------------------------------------
    %% T1w determinant
    % The determinant of the current matrix and old matrix should be the same,
    % otherwise this is suspicious of a left-right flip.

    Struct.ID = x.SUBJECTS{iSubject};

    PathOrientationResults = fullfile(x.dir.SUBJECTDIR,'xASL_qc_PrintOrientation_RigidRegT1.tsv');
    Struct.T1w_LR_flip_YesNo = uint8(xASL_im_DetermineFlip(x, iSubject, PathOrientationResults));
    % Whether left-right orientation has been flipped through registrations yes/no


    %% -----------------------------------------------------------------------------------------------
    %% WMH results
    PathLST = xASL_adm_GetFileList(x.D.TissueVolumeDir,['^WMH_LST_(LGA|LPA)_' Struct.ID '\.tsv$'],'FPList',[0 Inf]);

    if length(PathLST)>1
        warning('Too many LST volumetric results found, using first');
        bProcess = true;
    elseif length(PathLST)==1
        bProcess = true;        
    elseif isempty(PathLST)
        bProcess = false;
        if xASL_exist(x.P.Path_FLAIR,'file')
            warning('Didnt find any LST volumetric results');
            Struct.FLAIR_WMH_vol_mL = NaN;
            Struct.FLAIR_WMH_n = NaN;
            % WMH_vol_mL = White Matter Hyperintensity volume (mL)
            % WMH_n = Count of non-connected White Matter Hyperintensities
        end
    end
        
    if bProcess
        % Find WMH results in TSV
        [~, CellTSV] = xASL_bids_csv2tsvReadWrite(PathLST{1});
        for iC=1:size(CellTSV,1)
            if ~isempty(findstr(CellTSV{iC,1}, Struct.ID))
                Struct.FLAIR_WMH_vol_mL = CellTSV{iC,4};
                Struct.FLAIR_WMH_n = CellTSV{iC,5};
            end
        end    
    end


    %% -----------------------------------------------------------------------------------------------
    %% CAT12 QC output
    PathIQRresults = fullfile(x.D.TissueVolumeDir,['cat_' x.P.STRUCT '_' Struct.ID '.mat']);
    if exist(PathIQRresults,'file')
        curr = load(PathIQRresults);
        Struct.T1w_IQR_Perc = min(100,max(0,105 - curr.S.qualityratings.IQR*10));
        Struct.T1w_IQR_Perc = xASL_round(Struct.T1w_IQR_Perc,3);
    elseif xASL_exist(x.P.Path_T1,'file') && ~x.modules.structural.bSegmentSPM12
        warning('Didnt find any CAT12 volumetric results');
        Struct.T1w_IQR_Perc = NaN;
    end


    %% -----------------------------------------------------------------------------------------------
    %% CAT12 volumetric output
    PathCAT12Results = fullfile(x.D.TissueVolumeDir,['TissueVolume_' Struct.ID '.tsv']);
    if exist(PathCAT12Results,'file')
        [~, CellTSV] = xASL_bids_csv2tsvReadWrite(PathCAT12Results);

        for iC=1:size(CellTSV,1)
            if ~isempty(findstr(CellTSV{iC,1}, Struct.ID))
                Struct.T1w_GM_vol_mL = xASL_round(xASL_str2num(CellTSV{iC,2})*1000);
                Struct.T1w_WM_vol_mL = xASL_round(xASL_str2num(CellTSV{iC,3})*1000);
                Struct.T1w_CSF_vol_mL = xASL_round(xASL_str2num(CellTSV{iC,4})*1000,1); % bugfix SPM12 volume output

                Struct.T1w_ICV_vol_mL = xASL_round(Struct.T1w_GM_vol_mL+Struct.T1w_WM_vol_mL+Struct.T1w_CSF_vol_mL);
                Struct.T1w_GM_ICV_Ratio = xASL_round(Struct.T1w_GM_vol_mL/Struct.T1w_ICV_vol_mL,3);
                Struct.T1w_WM_ICV_Ratio = xASL_round(Struct.T1w_WM_vol_mL/Struct.T1w_ICV_vol_mL,3);
                Struct.T1w_CSF_ICV_Ratio = xASL_round(Struct.T1w_CSF_vol_mL/Struct.T1w_ICV_vol_mL,3);

            end
        end        
    else % fill in the blanks
        warning('Didnt find any CAT12 volumetric results');        
        Struct.T1w_GM_vol_mL = NaN;
        Struct.T1w_WM_vol_mL = NaN;
        Struct.T1w_CSF_vol_mL = NaN;
        Struct.T1w_ICV_vol_mL = NaN;
        Struct.T1w_GM_ICV_Ratio = NaN;
        Struct.T1w_WM_ICV_Ratio = NaN;
        Struct.T1w_CSF_ICV_Ratio = NaN;
    end

    %% -----------------------------------------------------------------------------------------------    
    %% Set struct fields to 4 decimals
    FieldNames  = fields(Struct);
    for iN=1:length(FieldNames)
        V = Struct.(FieldNames{iN});
        if isnumeric(V)
            Struct.(FieldNames{iN}) = xASL_round(V, 4);
        end
    end      

    x.Output.Structural = xASL_qc_FillFields(x.Output.Structural, Struct); % Fill fields

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