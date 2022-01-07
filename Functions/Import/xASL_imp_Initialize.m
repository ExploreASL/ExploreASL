function imPar = xASL_imp_Initialize(studyPath, imParPath)
%xASL_imp_Initialize Initialize the import.
%
% FORMAT: imPar = xASL_imp_Initialize(studyPath, imParPath)
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
% EXAMPLE:     imPar = xASL_imp_Initialize(studyPath, imParPath);
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% 1. Read study file

    % Initialize the imPar field
    [fpath, fname, fext] = fileparts(studyPath);

    % Load the imPar from the file
    if ~isempty(imParPath) && xASL_exist(imParPath,'file')==2
        % DCM2NII
        imPar = spm_jsonread(imParPath);
    else
        % NII2BIDS, DEFACE & BIDS2LEGACY
        imPar = struct;
    end

    %% 2. Specify paths
    if ~isfield(imPar, 'studyID') || isempty(imPar.studyID)
        imPar.studyID = [fname fext];
    end
    if ~isfield(imPar, 'TempRoot') || isempty(imPar.TempRoot) 
        imPar.TempRoot = fpath;
    end
    if ~isfield(imPar, 'RawRoot') || isempty(imPar.RawRoot)
        imPar.RawRoot = fpath;
    end
    if ~isfield(imPar, 'BidsRoot') || isempty(imPar.BidsRoot)
        imPar.BidsRoot = fpath;
    end

    %% 3. Finalize the directories
    imPar.RawRoot = fullfile(imPar.RawRoot, imPar.studyID, 'sourcedata');
    imPar.DerivativesRoot = fullfile(imPar.TempRoot, imPar.studyID, 'derivatives');
    imPar.TempRoot = fullfile(imPar.DerivativesRoot, 'ExploreASL', 'temp');
    imPar.LockRoot = fullfile(imPar.DerivativesRoot, 'ExploreASL', 'lock');
    imPar.BidsRoot = fullfile(imPar.BidsRoot, imPar.studyID, 'rawdata');

    %% 4. Specify the tokens
    
    % Specify the imPar struct
    imPar = xASL_imp_InitializeBuildImPar(imPar);
    
    % Warn the user if token aliases are invalid according to our pipeline/BIDS
    xASL_imp_InitializeCheckTokens(imPar);

    %% 5. Specify the additional details of the conversion
    if ~isfield(imPar, 'bVerbose') || isempty(imPar.bVerbose)
        imPar.bVerbose = true;
    end
    if ~isfield(imPar, 'bOverwrite') || isempty(imPar.bOverwrite)
        imPar.bOverwrite  = false; % NB, the summary file will be recreated anyway and dicom conversion in temp is always done, even if dest. exists
    end
    if ~isfield(imPar ,'visitNames') || isempty(imPar.visitNames)
        imPar.visitNames = {};
    end
    if ~isfield(imPar,'sessionNames') || isempty(imPar.sessionNames)
        imPar.sessionNames = {};
    end
    if ~isfield(imPar, 'sessionNames') || isempty(imPar.sessionNames)
        imPar.sessionNames = {};
    end
    
    % Create the empty derivatives directory for the general import (we need it for lock and temp files)
    xASL_adm_CreateDir(imPar.DerivativesRoot);
    
    % Create the lock files directory (especially for the import lock files)
    xASL_adm_CreateDir(imPar.LockRoot);

end

%% Specify the imPar struct
function imPar = xASL_imp_InitializeBuildImPar(imPar)

    if ~isfield(imPar, 'folderHierarchy')
        imPar.folderHierarchy = {}; % must define this per study; use a cell array of regular expressions. One cell per directory level.
    end
    if ~isfield(imPar, 'tokenOrdering')
        imPar.tokenOrdering = []; % must match imPar.folderHierarchy: 1==subject, 2=visit, 3==session, 4==scan (if visit or session are omitted, they will be skipped)
    else
        imPar.tokenOrdering = imPar.tokenOrdering(:)';
    end
    if ~isfield(imPar, 'tokenScanAliases')
        imPar.tokenScanAliases = [];
    else
        if (size(imPar.tokenScanAliases,2) > 2) || (size(imPar.tokenScanAliases,2) == 1)
            tokenScanAliasesOld = imPar.tokenScanAliases;
			imPar = rmfield(imPar, 'tokenScanAliases');
            imPar.tokenScanAliases(:,1) = tokenScanAliasesOld(1:2:end);
            imPar.tokenScanAliases(:,2) = tokenScanAliasesOld(2:2:end);
        end
    end
    if ~isfield(imPar, 'tokenVisitAliases')
        imPar.tokenVisitAliases = [];
    else
        if (size(imPar.tokenVisitAliases,2) > 2) || (size(imPar.tokenVisitAliases,2) == 1)
            tokenVisitAliasesOld = imPar.tokenVisitAliases;
			imPar = rmfield(imPar, 'tokenVisitAliases');
            imPar.tokenVisitAliases(:,1) = tokenVisitAliasesOld(1:2:end);
            imPar.tokenVisitAliases(:,2) = tokenVisitAliasesOld(2:2:end);
        end
    end
    if ~isfield(imPar, 'tokenSessionAliases')
        imPar.tokenSessionAliases = [];
    else
        if (size(imPar.tokenSessionAliases,2) > 2) || (size(imPar.tokenSessionAliases,2) == 1)
            tokenSessionAliasesOld = imPar.tokenSessionAliases;
			imPar = rmfield(imPar, 'tokenSessionAliases');
            imPar.tokenSessionAliases(:,1) = tokenSessionAliasesOld(1:2:end);
            imPar.tokenSessionAliases(:,2) = tokenSessionAliasesOld(2:2:end);
        end
    end
    if ~isfield(imPar,'bMatchDirectories')
        imPar.bMatchDirectories  = false;
    end

end


%% Warn the user if token aliases are invalid according to our pipeline/BIDS
function xASL_imp_InitializeCheckTokens(imPar)

    % Check tokenScanAliases
    if isfield(imPar, 'tokenScanAliases')
        if size(imPar.tokenScanAliases,2)==2
            for iToken = 1:size(imPar.tokenScanAliases,1)
                currentToken = imPar.tokenScanAliases{iToken,2};
                % Warn if the token contains ASL but is not ASL4D
                if ~isempty(regexpi(currentToken,'ASL', 'once')) && isempty(regexp(currentToken, 'ASL4D', 'once'))
                    warning('Please use ASL4D for your ASL tokenScanAlias instead...');
                end
                % Warn if the token contains T1 but is not T1w
                if ~isempty(regexpi(currentToken, 'T1', 'once')) && ~strcmp(currentToken, 'T1w')
                    warning('Please use T1w for your T1 tokenScanAlias instead...');
                end
            end
        end
    end

end

