function [x] = xASL_init_BIDS2Legacy(x)
%xASL_init_BIDS2Legacy Prepare relevant fields in xASL x struct for BIDS to Legacy conversion
%
% FORMAT: [x] = xASL_init_BIDS2Legacy(x)
% 
% INPUT:
%   x          - ExploreASL x structure (STRUCT, REQUIRED)
%
% OUTPUT:
%   x          - ExploreASL x structure
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Checks the necessary input directories and loads the BIDS directory structure that should be done only once for all subjects and not
%                 for each subject separately
%
% 1. Check basic directories
% 2. Start with checking dataset_description.json & rawdata
% - 1. The input is dataset_description.json in the rawdata folder
% - 2. The input is dataPar.json or sourceStructure.json - have to look for a rawdata folder
% 3. Check if a dataPar is provided, otherwise use the defaults
% 4. Load the BIDS structure of all subjects
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        [x] = xASL_init_BIDS2Legacy(x);
% __________________________________
% Copyright (c) 2015-2023 ExploreASL

% Admin
if nargin < 1 || isempty(x)
	error('x-struct is a required input');
end

% Update paths
x = xASL_imp_UpdateDatasetRoot(x);


%% 1. Check basic directories
if ~isfield(x,'dir')
	error('Missing directories field...');
end

if ~isfield(x.dir,'dataset_description')
	error('Missing dataset_description field...');
end
if isempty(x.dir.dataset_description)
	error('Empty dataset_description path...');
end

if isempty(x.dir.DatasetRoot)
	error('x.dir.DatasetRoot is a required parameter.');
end

% Verify that the rawdata subfolder exists
if ~exist(fullfile(x.dir.DatasetRoot,'rawdata'),'dir')
    warning('Invalid folder selected, not containing rawdata folder');
    return;
end

% Check derivatives of ExploreASL
if ~isfield(x.dir,'xASLDerivatives')
    error('Missing xASL derivatives field...');
end


%% 2. Start with checking dataset_description.json & rawdata
% Either the rawdata/dataset_description.json was provided
% Or the dataPar.json was provided and we have to search was dataset_description.json

[rawdataFolder, fileNameDatasetDescription, fileExtension] = fileparts(x.dir.dataset_description);
if strcmpi([fileNameDatasetDescription fileExtension], 'dataset_description.json')
	% 2.1 The input is dataset_description.json in the rawdata folder
	[~, rawdataFolderName] = fileparts(rawdataFolder);
	if ~strcmp(rawdataFolderName, 'rawdata')
		error('Invalid folder in which dataset_description.json was found, should be /rawdata');
	end
else
	% 2.2 The input is dataPar.json or sourceStructure.json - have to look for a rawdata folder
	if ~exist(x.dir.RawData,'dir')
		error('Path rawdata does not exist...');
	elseif ~exist(fullfile(x.dir.RawData,'dataset_description.json'),'file')
		error('File dataset_description.json is not found...');
	end
end


%% 3.Check if a dataPar is provided, otherwise use the defaults

% Get dataPar
fListDataPar = xASL_adm_GetFileList(x.dir.DatasetRoot,'(?i)(^dataPar.*\.json$)', 'FPList', [], 0);
if length(fListDataPar) < 1
	% Fill the dataPars with default parameters
	fprintf('There is no dataPar.json file in the study root directory. Default settings will be used...\n');
else
	if length(fListDataPar) > 1
		warning('Multiple dataPar*.jsons exist. Using the first: %s\n',fListDataPar{1});
	end
	% Fill the dataPars with the provided parameters
	x.dataPar = xASL_io_ReadJson(fListDataPar{1});
end


%% 4. Load the BIDS structure of all subjects
x.modules.bids2legacy.BIDS = bids.layout(fullfile(x.dir.DatasetRoot,'rawdata'));


%% 5. Create the derivatives directory
if exist(x.dir.xASLDerivatives, 'dir')
    fprintf('%s\n', [x.dir.xASLDerivatives ' already exists']);
    fprintf('%s\n', 'Note that all pre-existing derivative subject folders will be overwritten,');
    fprintf('%s\n', 'unless BIDS2Legacy lock files already exist for a subject');
else
    xASL_adm_CreateDir(x.dir.xASLDerivatives);
end


%% 6. Load BIDS configuration for file renaming
x.modules.bids2legacy.bidsPar = xASL_bids_Config;


    
% This is the line needed by xASL_init_Iteration for BIDS2Legacy
x.D.ROOT = x.dir.DatasetRoot;

% SESSIONS DUMMY
x.SESSIONS = {''};


end