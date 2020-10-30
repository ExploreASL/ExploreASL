function [identical,results] = xASL_bids_compareStructures(pathDatasetA,pathDatasetB)
%xASL_bids_compareStructures Function that compares two BIDS folders with several subfolders and studies and prints the differences.
%
% FORMAT: [identical,results] = xASL_bids_compareStructures(pathDatasetA,pathDatasetB);
%
% INPUT:
%        pathDatasetA       - path to first BIDS structure (REQUIRED)
%        pathDatasetB       - path to second BIDS structure (REQUIRED)
%
% OUTPUT:
%        identical          - Returns 1 if both folder structures are identical and 0 if not
%        results            - structure containing (possible) differences of both folder structures
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      Function that compares two BIDS folders with several subfolders and studies and prints the differences.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          pathDatasetA = '...\bids-examples\eeg_rest_fmri';
%                   pathDatasetB = '...\bids-examples\eeg_rest_fmri_exact_copy'
%                   [identical,results] = xASL_bids_compareStructures(pathDatasetA,pathDatasetB);
%
% REFERENCES:       ...
% __________________________________
% Copyright @ 2015-2020 ExploreASL


%% Input Check

% Check if both root folders are valid char arrays or strings
if ~(ischar(pathDatasetA) || isstring(pathDatasetA))
    error('The path of structure A is neither a char array not a string...');
end
if ~(ischar(pathDatasetB) || isstring(pathDatasetB))
    error('The path of structure A is neither a char array not a string...');
end

% Check if both root folders exists
if ~(xASL_exist(pathDatasetA)==7)
    error('The root folder of structure A does not exist...');
end
if ~(xASL_exist(pathDatasetB)==7)
    error('The root folder of structure B does not exist...');
end


%% Defaults

% Set identical to true (will be set to false as soon as a difference is found)
identical = true;

% Initialize results structure
results = struct;

%% Check first level (expected files and folders: participants, dataset_description, sub-xx, code, sourcedata, derivatives etc.)

% Get dataset names
[~,datasetA,~] = fileparts(pathDatasetA);
[~,datasetB,~] = fileparts(pathDatasetB);

results.(datasetA) = struct;
results.(datasetB) = struct;

% Get all files and folders of dataset A (level 1)
datasetA_files_level1 = xASL_adm_GetFileList(pathDatasetA,[],false);
datasetA_folders_level1 = xASL_adm_GetFileList(pathDatasetA,[],false,[],true);

% Get all files and folders of dataset B (level 1)
datasetB_files_level1 = xASL_adm_GetFileList(pathDatasetB,[],false);
datasetB_folders_level1 = xASL_adm_GetFileList(pathDatasetB,[],false,[],true);

% Check if all files of A are in B and the other way round
if ~isempty(setdiff(datasetA_files_level1,datasetB_files_level1)) || ~isempty(setdiff(datasetB_files_level1,datasetA_files_level1))
    % Set identical to false
    identical = false;
    fprintf(repmat('=',100,1));
    % Get missing files of A that are not in B
    fprintf('\nLEVEL 1: %s\n',datasetA);
    results.(datasetB).missingFiles = setdiff(datasetA_files_level1,datasetB_files_level1);
    if ~isempty(results.(datasetB).missingFiles)
        % Iterate over files
        for it=1:length(results.(datasetB).missingFiles)
            fprintf('FILE:    %s\n',results.(datasetB).missingFiles{it,1});
        end
    end
    % Get missing files of B that are not in A
    fprintf('LEVEL 1: %s\n',datasetB);
    results.(datasetA).missingFiles = setdiff(datasetB_files_level1,datasetA_files_level1);
    if ~isempty(results.(datasetA).missingFiles)
        % Iterate over files
        for it=1:length(results.(datasetA).missingFiles)
            fprintf('FILE:    %s\n',results.(datasetA).missingFiles{it,1});
        end
    end
end

% Check if all folders of A are in B and the other way round
if ~isempty(setdiff(datasetA_folders_level1,datasetB_folders_level1)) || ~isempty(setdiff(datasetB_folders_level1,datasetA_folders_level1))
    % Set identical to false
    identical = false;
    fprintf(repmat('=',100,1));
    % Get missing folders of A that are not in B
    fprintf('\nLEVEL 1: %s\n',datasetA);
    results.(datasetB).missingFolders = setdiff(datasetA_folders_level1,datasetB_folders_level1);
    if ~isempty(results.(datasetB).missingFolders)
        % Iterate over files
        for it=1:length(results.(datasetB).missingFolders)
            fprintf('FILE:    %s\n',results.(datasetB).missingFolders{it,1});
        end
    end
    % Get missing folders of B that are not in A
    fprintf('LEVEL 1: %s\n',datasetB);
    results.(datasetA).missingFolders = setdiff(datasetB_folders_level1,datasetA_folders_level1);
    if ~isempty(results.(datasetA).missingFolders)
        % Iterate over files
        for it=1:length(results.(datasetA).missingFolders)
            fprintf('FILE:    %s\n',results.(datasetA).missingFolders{it,1});
        end
    end
end

