function [jsonLocal, bJsonLocalM0isFile] = xASL_imp_NII2BIDS_Subject_DefineM0Type(studyPar, bidsPar, jsonLocal, inSessionPath, subjectLabel, sessionLabel)
%xASL_imp_NII2BIDS_Subject_DefineM0Type Define M0 Type.
%
% FORMAT: [jsonLocal, bJsonLocalM0isFile] = xASL_imp_NII2BIDS_Subject_DefineM0Type(studyPar, bidsPar, jsonLocal, inSessionPath, subjectLabel, sessionLabel)
% 
% INPUT:
%   studyPar       - JSON file with the BIDS parameters relevant for the whole study (STRUCT, REQUIRED)
%   bidsPar        - Output of xASL_imp_Config (STRUCT, REQUIRED)
%   jsonLocal      - jsonLocal struct (STRUCT, REQUIRED)
%   inSessionPath  - inSession path (CHAR ARRAY, PATH, REQUIRED)
%   subjectLabel   - subject label (CHAR ARRAY, REQUIRED)
%   sessionLabel   - session label (CHAR ARRAY, REQUIRED)
%
% OUTPUT:
%   jsonLocal      - jsonLocal struct (STRUCT, REQUIRED)
%   bJsonLocalM0isFile    - JSON local M0 is file (BOOLEAN)
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Define M0 Type.
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     [jsonLocal, bJsonLocalM0isFile] = xASL_imp_NII2BIDS_Subject_DefineM0Type(studyPar, bidsPar, jsonLocal, inSessionPath, subjectLabel, sessionLabel);
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% Type of an M0 image
    bJsonLocalM0isFile = 0;
    if ~isfield(studyPar,'M0') || isempty(studyPar.M0) || strcmpi(studyPar.M0,'separate_scan')
        if isfield(studyPar,'M0PositionInASL4D') && (max(studyPar.M0PositionInASL4D(:))>0)
            jsonLocal.M0 = true;
            jsonLocal.M0Type = bidsPar.strM0Included;
        elseif xASL_exist(fullfile(inSessionPath,'M0.nii'))
            jsonLocal.M0 = fullfile(bidsPar.strPerfusion,['sub-' subjectLabel sessionLabel]);
            jsonLocal.M0Type = bidsPar.strM0Separate;
            bJsonLocalM0isFile = 1;
        else
            if ~isempty(strfind(jsonLocal.ASLContext,bidsPar.strM0scan))
                jsonLocal.M0 = true;
                jsonLocal.M0Type = bidsPar.strM0Included;
            else
                jsonLocal.M0 = false;
                jsonLocal.M0Type = bidsPar.strM0Absent;
            end
        end
    else
        if strcmpi(studyPar.M0,'UseControlAsM0')
            jsonLocal.M0 = bidsPar.strM0Absent;
        else
            if strcmpi(studyPar.M0,'no_background_suppression')
                jsonLocal.M0 = bidsPar.strM0Absent;
            else
                jsonLocal.M0 = studyPar.M0;
                if isnumeric(studyPar.M0)
                    jsonLocal.M0Type = bidsPar.strM0Estimate;
                    jsonLocal.M0Estimate = studyPar.M0;
                elseif xASL_exist(fullfile(inSessionPath,'M0.nii'))
                    jsonLocal.M0Type = bidsPar.strM0Separate;
                elseif ~isfield(jsonLocal, 'ASLContext')
                    warning('jsonLocal.ASLContext missing, this may crash');
                elseif ~isempty(strfind(jsonLocal.ASLContext,bidsPar.strM0scan))
                    jsonLocal.M0Type = bidsPar.strM0Included;
                else
                    jsonLocal.M0Type = bidsPar.strM0Absent;
                end
            end
        end
    end
    
end


