function [x] = xASL_init_Loaddataset_description(x)
%xASL_init_Loaddataset_description Load data parameter file
%
% FORMAT: [x] = xASL_init_Loaddataset_description(x)
%
% INPUT:
%   x       - ExploreASL x structure (STRUCT, REQUIRED)
%
% OUTPUT:
%   x       - ExploreASL x structure (STRUCT)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Load dataset_decription.json file (which contains BIDS parameters about the dataset in rawdata
%
% EXAMPLE:     This is part of the initialization workflow. Check out the usage there.
%
% This function performs the following parts:
% 1. Find the dataset_description.json location
% 2. Select dataset_description location
% 3. Issue warnings if we have multiple dataset_descriptives.json files at the selected location
% 2. Move dataset_description.json to derivatives
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:  n/a
% __________________________________
% Copyright (c) 2015-2024 ExploreASL


%% 1. Find the dataset_description.json location
% Any of these bUse booleans becomes true if we find a dataset_description.json file on the respective location

bUseRoot = false;
bUseRawdata = false;
bUseDerivatives = false;

% Check dataset_description inside the root folder
listRoot = xASL_adm_GetFileList(x.dir.DatasetRoot, '(?i)(^dataset_description.*\.json$)', 'FPList', [], 0);
if ~isempty(listRoot)
    warning([xASL_num2str(length(listRoot)) ' dataset_description.json (or similar) file(s) found in ' x.dir.DatasetRoot ', will try this but this is not the appropriate location']);
    bUseRoot = true;
end   

% Check dataset_description inside the rawdata folder /rawdata
listRawdata = xASL_adm_GetFileList(x.dir.RawData, '(?i)(^dataset_description.*\.json$)', 'FPList', [], 0);
if ~isempty(listRawdata)
    bUseRawdata = true;
else
    warning('No valid dataset_description.json found within the rawdata directory...');
    fprintf('%s\n', 'Ensure that your data is BIDS compatible by running the BIDS validator');
end

% Check dataset_description inside the processing folder /derivatives/ExploreASL
fListLegacy = xASL_adm_GetFileList(x.dir.xASLDerivatives, '(?i)(^dataset_description.*\.json$)', 'FPList', [], 0);
if ~isempty(fListLegacy)
    bUseDerivatives = true;
end


%% 2. Select dataset_description location
% Here we check which locations to use based on the bUse booleans
% We prioritize according to:
% 1) /derivatives/ExploreASL (not replacing existing derivatives/ExploreASL files)
% 2) /rawdata (this is where the file should be according to BIDS)
% 3) /root

% Prioritization
if bUseDerivatives && bUseRawdata
    fprintf('%s\n', 'dataset_description.json (or similar) file(s) both in /derivatives/ExploreASL and /rawdata folder, ignoring the last');
    bUseRawdata = false;
end
if bUseRawdata && bUseRoot
    fprintf('%s\n', 'dataset_description.json (or similar) file(s) both in /rawdata and / root folder, ignoring the last');
    bUseRoot = false;
end

% Select the folder
if bUseDerivatives
    fFolder = x.dir.xASLDerivatives;
    listDataSetDescript = fListLegacy;
    fprintf('%s\n', ['dataset_descriptives.json found in ' fFolder]);
elseif bUseRawdata
    fFolder = x.dir.RawData;
    listDataSetDescript = listRawdata;
    fprintf('%s\n', ['dataset_descriptives.json found in ' fFolder]);
elseif bUseRoot
    fFolder = x.dir.DatasetRoot;
    listDataSetDescript = listRoot;
    fprintf('%s\n', ['dataset_descriptives.json found in ' fFolder]);
else
    % Issue a warning if no files are found, and abort this function
    warning('No dataset_descriptives.json found! Check if you data are BIDS compatible using the BIDS validator');
    return;
end


%% 3. Issue warnings if we have multiple dataset_descriptives.json files at the selected location
if length(listDataSetDescript)>1
    fprintf('Warning: multiple dataset_descriptives.json files found, using the first\n');
end


%% 4. Move dataset_description.json to derivatives
destinationFile = fullfile(x.dir.xASLDerivatives, 'dataset_description.json');
if ~xASL_exist(destinationFile)
    % We prioritize pre-existing files in /derivatives/ExploreASL
    xASL_Copy(listDataSetDescript{1}, destinationFile);
end


end