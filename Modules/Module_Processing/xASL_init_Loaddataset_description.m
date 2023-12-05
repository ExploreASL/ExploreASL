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
% 1. Admin checks
% 2. Find the dataset_description.json location
% 3. Load pre-existing dataset_descriptives.json
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:  n/a
% __________________________________
% Copyright (c) 2015-2024 ExploreASL


%% 1. Admin checks
if ~isfield(x.dir,'dataset_description')
	error('Missing dataset_description field...');
end
if isempty(x.dir.dataset_description)
	error('Empty dataset_description path...');
end


%% 2. Find the dataset_description.json location

bUseRoot = false;
bUseRawdata = false;
bUseDerivatives = false;

% Check dataPar inside the root folder
listRoot = xASL_adm_GetFileList(x.dir.DatasetRoot, '(?i)(^dataset_description.*\.json$)', 'FPList', [], 0);
if ~isempty(listRoot)
    warning([xASL_num2str(length(listRoot)) ' dataset_description.json (or similar) file(s) found in ' x.dir.DatasetRoot ', will try this but this is not the appropriate location']);
    bUseRoot = true;
end   

% Check dataPar inside the rawdata folder /rawdata
listRawdata = xASL_adm_GetFileList(x.dir.RawData, '(?i)(^dataset_description.*\.json$)', 'FPList', [], 0);
if ~isempty(listRawdata)
    bUseRawdata = true;
else
    warning('No valid dataset_description.json found within the rawdata directory...');
    fprintf('%s\n', 'Ensure that your data is BIDS compatible by running the BIDS validator');
end

% Check dataPar inside the processing folder /derivatives/ExploreASL
fListLegacy = xASL_adm_GetFileList(x.dir.xASLDerivatives, '(?i)(^dataset_description.*\.json$)', 'FPList', [], 0);
if ~isempty(fListLegacy)
    bUseDerivatives = true;
end

% Choose folder & filelist
if bUseDerivatives && bUseRawdata
    % We prioritize existing dataset_description.json in the derivatives folder
    fprintf('%s\n', 'dataset_description.json (or similar) file(s) both in /derivatives/ExploreASL and /rawdata folder, ignoring the last');
    bUseRawdata = false;
end

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
    listDataSetDescript = {};
end

destinationFile = fullfile(x.dir.xASLDerivatives, 'dataset_description.json');
if ~xASL_exist(destinationFile)
    xASL_Copy(listDataSetDescript{1}, destinationFile);
end


%% 3. Load pre-existing dataset_descriptives.json
if length(listDataSetDescript)>1
    fprintf('Warning: multiple dataset_descriptives.json files found\n');
elseif isempty(listDataSetDescript)
    warning('No dataset_descriptives.json provided');
end


end