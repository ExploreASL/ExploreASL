function [x] = xASL_imp_BIDS2Legacy(x)
%xASL_imp_BIDS2Legacy BIDS2LEGACY conversion script which calls xASL_bids_BIDS2Legacy.
%
% FORMAT: [x] = xASL_imp_BIDS2Legacy(x);
%
% INPUT:
%   x             - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x             - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    BIDS to Legacy conversion script which calls xASL_bids_BIDS2Legacy.
%
% 1. Input check
% 2. Start with checking dataset_description.json & rawdata
% - 1. The input is dataset_description.json in the rawdata folder
% - 2. The input is dataPar.json or sourceStructure.json - have to look for a rawdata folder
% 3. Run the legacy conversion: Check if a dataPar is provided, otherwise use the defaults
% 4. Overwrite DatasetRoot
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% 1. Input check
    if ~isfield(x,'dir')
        error('Missing directories field...');
    end
    if ~isfield(x.dir,'dataset_description')
        error('Missing dataset_description field...');
    end
    if isempty(x.dir.dataset_description)
        error('Empty dataset_description path...');
    end


    %% 2. Start with checking dataset_description.json & rawdata
    % Either the rawdata/dataset_description.json was provided
	% Or the dataPar.json was provided and we have to search was dataset_description.json
    
    [rawdataFolder, fileNameDatasetDescription, fileExtension] = fileparts(x.dir.dataset_description);
    if strcmpi([fileNameDatasetDescription fileExtension], 'dataset_description.json')
        %% 2.1 The input is dataset_description.json in the rawdata folder
        [~, rawdataFolderName] = fileparts(rawdataFolder);
        if ~strcmp(rawdataFolderName, 'rawdata')
            error('Invalid folder in which dataset_description.json was found, should be /rawdata');
        end
        localDatasetRoot = x.dir.DatasetRoot;
    else
        %% 2.2 The input is dataPar.json or sourceStructure.json - have to look for a rawdata folder
        pathRawData = fullfile(x.dir.DatasetRoot,'rawdata');
        % Check the if the correct BIDS structure is in place
        if ~exist(pathRawData,'dir')
            error('Path rawdata does not exist');
        elseif ~exist(fullfile(pathRawData,'dataset_description.json'),'file')
            error('File dataset_description.json is not found');
        else
            localDatasetRoot = x.dir.DatasetRoot;
        end
    end

	%% 3. Run the legacy conversion: Check if a dataPar is provided, otherwise use the defaults
	fListDataPar = xASL_adm_GetFileList(localDatasetRoot,'(?i)(^dataPar.*\.json$)', 'FPList', [], 0);
	if length(fListDataPar) < 1
		fprintf('There is no dataPar.json file in the study root directory. Default settings will be used...\n');
		% Fill the dataPars with default parameters
		dataPar = xASL_bids_BIDS2Legacy(localDatasetRoot, x, 1, []);
	else
		% Fill the dataPars with the provided parameters
		dataPar = spm_jsonread(fListDataPar{1});
		dataPar = xASL_bids_BIDS2Legacy(localDatasetRoot, x, 1, dataPar);
	end

	%% 4. Overwrite DatasetRoot
    fieldsDataPar = fieldnames(dataPar.x);
    % Add fields that are in dataPar.x but missing in x
    for iField = 1:numel(fieldsDataPar)
        if ~isfield(x,fieldsDataPar{iField,1}) && ~strcmp('dir',fieldsDataPar{iField,1})
            x.(fieldsDataPar{iField,1}) = dataPar.x.(fieldsDataPar{iField,1});
        end
    end
    % Update dataPar path
	x.dir.dataPar = dataPar.x.dir.dataPar;
    
end


