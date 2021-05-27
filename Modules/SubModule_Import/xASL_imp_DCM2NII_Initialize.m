function imPar = xASL_imp_DCM2NII_Initialize(studyPath, imParPath)
%xASL_imp_DCM2NII_Initialize Initialize DCM2NII.
%
% FORMAT: imPar = xASL_imp_DCM2NII_Initialize(studyPath, imParPath)
% 
% INPUT:
%   studyPath  - Path to the study directory containing the 'sourcedata' directory with the DICOM files (REQUIRED, CHAR ARRAY)
%   imParPath  - Path to the JSON file with structure with import parameters (REQUIRED, CHAR ARRAY)
%
% OUTPUT:
%   imPar      - JSON file with structure with import parameters (REQUIRED, STRUCT)
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Initialize DCM2NII.
%
% 1. Read study file
% 2. Specify paths
% 3. Finalize the directories
% 4. Specify the tokens
% 5. Specify the additional details of the conversion
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     imPar = xASL_imp_DCM2NII_Initialize(studyPath, imParPath);
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% 1. Read study file

    % Initialize the imPar field
    [fpath, fname, fext] = fileparts(studyPath);

    % Load the imPar from the file
    imPar = spm_jsonread(imParPath);

    %% 2. Specify paths
    if ~isfield(imPar,'studyID') || isempty(imPar.studyID)
        imPar.studyID = [fname fext];
    end
    if ~isfield(imPar,'TempRoot') || isempty(imPar.TempRoot) 
        imPar.TempRoot = fpath;
    end
    if ~isfield(imPar,'RawRoot') || isempty(imPar.RawRoot)
        imPar.RawRoot = fpath;
    end
    if ~isfield(imPar,'BidsRoot') || isempty(imPar.BidsRoot)
        imPar.BidsRoot = fpath;
    end

    %% 3. Finalize the directories
    imPar.RawRoot = fullfile(imPar.RawRoot,imPar.studyID, 'sourcedata'); % default name
    imPar.TempRoot = fullfile(imPar.TempRoot,imPar.studyID,'temp');
    imPar.BidsRoot = fullfile(imPar.BidsRoot,imPar.studyID,'rawdata');

    %% 4. Specify the tokens
    if ~isfield(imPar,'folderHierarchy')
        imPar.folderHierarchy = {}; % must define this per study; use a cell array of regular expressions. One cell per directory level.
    end
    if ~isfield(imPar,'tokenOrdering')
        imPar.tokenOrdering = []; % must match imPar.folderHierarchy: 1==subject, 2=visit, 3==session, 4==scan (if visit or session are omitted, they will be skipped)
    else
        imPar.tokenOrdering = imPar.tokenOrdering(:)';
    end
    if ~isfield(imPar,'tokenScanAliases')
        imPar.tokenScanAliases = [];
    else
        if (size(imPar.tokenScanAliases,2) > 2) || (size(imPar.tokenScanAliases,2) == 1)
            tokenScanAliasesOld = imPar.tokenScanAliases;
			imPar = rmfield(imPar,'tokenScanAliases');
            imPar.tokenScanAliases(:,1) = tokenScanAliasesOld(1:2:end);
            imPar.tokenScanAliases(:,2) = tokenScanAliasesOld(2:2:end);
        end
    end
    if ~isfield(imPar,'tokenVisitAliases')
        imPar.tokenVisitAliases = [];
    else
        if (size(imPar.tokenVisitAliases,2) > 2) || (size(imPar.tokenVisitAliases,2) == 1)
            tokenVisitAliasesOld = imPar.tokenVisitAliases;
			imPar = rmfield(imPar,'tokenVisitAliases');
            imPar.tokenVisitAliases(:,1) = tokenVisitAliasesOld(1:2:end);
            imPar.tokenVisitAliases(:,2) = tokenVisitAliasesOld(2:2:end);
        end
    end
    if ~isfield(imPar,'tokenSessionAliases')
        imPar.tokenSessionAliases = [];
    else
        if (size(imPar.tokenSessionAliases,2) > 2) || (size(imPar.tokenSessionAliases,2) == 1)
            tokenSessionAliasesOld = imPar.tokenSessionAliases;
			imPar = rmfield(imPar,'tokenSessionAliases');
            imPar.tokenSessionAliases(:,1) = tokenSessionAliasesOld(1:2:end);
            imPar.tokenSessionAliases(:,2) = tokenSessionAliasesOld(2:2:end);
        end
    end
    if ~isfield(imPar,'bMatchDirectories')
        imPar.bMatchDirectories  = false;
    end

    %% 5. Specify the additional details of the conversion
    if ~isfield(imPar,'bVerbose') || isempty(imPar.bVerbose)
        imPar.bVerbose = true;
    end
    if ~isfield(imPar,'bOverwrite') || isempty(imPar.bOverwrite)
        imPar.bOverwrite  = false; % NB, the summary file will be recreated anyway and dicom conversion in temp is always done, even if dest. exists
    end
    if ~isfield(imPar,'visitNames') || isempty(imPar.visitNames)
        imPar.visitNames = {};
    end
    if ~isfield(imPar,'nMaxVisits') || isempty(imPar.nMaxVisits)
        imPar.nMaxVisits = 0;
    end
    if ~isfield(imPar,'sessionNames') || isempty(imPar.sessionNames)
        imPar.sessionNames = {};
    end
    if ~isfield(imPar,'nMaxSessions') || isempty(imPar.nMaxSessions)
        imPar.nMaxSessions = 0;
    end

end




