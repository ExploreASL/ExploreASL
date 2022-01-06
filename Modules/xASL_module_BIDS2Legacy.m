function [result, x] = xASL_module_BIDS2Legacy(x)
%xASL_module_BIDS2Legacy BIDS2LEGACY conversion script which calls xASL_bids_BIDS2Legacy.
%
% FORMAT: [x] = xASL_module_BIDS2Legacy(x);
%
% INPUT:
%   x             - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x             - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    BIDS to Legacy conversion script which calls xASL_module_BIDS2Legacy.
%
% 1. Input check
% 2. Start with checking dataset_description.json & rawdata
% - 1. The input is dataset_description.json in the rawdata folder
% - 2. The input is dataPar.json or sourceStructure.json - have to look for a rawdata folder
% 3. Run the legacy conversion: Check if a dataPar is provided, otherwise use the defaults
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
%
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% 1. Input check
    
    % Make sure that logging is still active
    if isfield(x.dir,'diaryFile')
        diary(x.dir.diaryFile);
    end
    
    % Start Mutex
    x = xASL_init_InitializeMutex(x, 'BIDS2Legacy');
    
    % Define lock states
    StateName{1} = '010_BIDS2LEGACY';
    StateName{2} = '020_CleanUp';
    
    % Default for result
    result = true;
    
    % Print feedback
    xASL_adm_BreakString('BIDS to ExploreASL LEGACY CONVERSION');
    
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
        % 2.1 The input is dataset_description.json in the rawdata folder
        [~, rawdataFolderName] = fileparts(rawdataFolder);
        if ~strcmp(rawdataFolderName, 'rawdata')
            error('Invalid folder in which dataset_description.json was found, should be /rawdata');
        end
    else
        % 2.2 The input is dataPar.json or sourceStructure.json - have to look for a rawdata folder
        pathRawData = fullfile(x.dir.DatasetRoot,'rawdata');
        % Check the if the correct BIDS structure is in place
        if ~exist(pathRawData,'dir')
            error('Path rawdata does not exist');
        elseif ~exist(fullfile(pathRawData,'dataset_description.json'),'file')
            error('File dataset_description.json is not found');
        end
    end

	%% 3. Run the legacy conversion: Check if a dataPar is provided, otherwise use the defaults
    
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
        x.dataPar = spm_jsonread(fListDataPar{1});
    end
    % Run legacy conversion
    iState = 1;
    if ~x.mutex.HasState(StateName{1})
        x = xASL_wrp_BIDS2Legacy(x, 1);
        x.mutex.AddState(StateName{iState});
    elseif x.mutex.HasState(StateName{1})
        fprintf('BIDS2Legacy was run before...   \n');
    end
    
    %% Clean-Up
    x = xASL_imp_CompleteBIDS2Legacy(x);
    x.mutex.AddState('999_ready');
    x.mutex.Unlock();
    
    
end


