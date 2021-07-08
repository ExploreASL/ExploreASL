function [x] = xASL_init_checkDatasetRoot(x)
%xASL_init_checkDatasetRoot Check the ExploreASL parameter "DatasetRoot"
%
% FORMAT: 
%   [x] = xASL_init_checkDatasetRoot(x)
%
% INPUT:
%   x             - ExploreASL x structure (REQUIRED, STRUCT)
%
% OUTPUT:
%   x             - ExploreASL x structure
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Check the ExploreASL parameter "DatasetRoot".
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     n/a
%
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% Check the ExploreASL parameter "DatasetRoot"

    % Default
    if ~isfield(x.opts,'bLoadData')
        x.opts.bLoadData = false;
    end

    % Check if the DatasetRoot is a directory (NEW - ASL BIDS)
    x.opts.dataParType = 'unknown'; % default
    % Create directory field if it doesn't exist already
    if ~isfield(x, 'dir')
        x.dir = struct;
    end

    % Check if user correctly inserted a dataset root directory
    if exist(x.opts.DatasetRoot,'dir')
        [x] = xASL_init_DetermineRequiredPaths(x);
    elseif exist(x.opts.DatasetRoot,'file')
        % Temporary functionality, this will lead to an error starting v2.0.0
        [x] = xASL_init_checkDatasetRoot_invalid_starting_2_0(x);
    else
        if ~isempty(x.opts.DatasetRoot) % The user inserted a directory or file which does not exist
            warning('Dataset root directory does not exist...');
        end
        if x.opts.bProcessData || x.opts.bImportData
            if ~isdeployed
                x.opts.DatasetRoot = input('Please insert the path to your study directory: ');
            else
                error('Study directory does not exist...');
            end
            % Immediately check the input
            if ~exist(x.opts.DatasetRoot, 'dir')
                warning('This study directory does not exist, ExploreASL will only be initialized...');
                x.opts.bProcessData = 0;
                x.opts.bImportData = 0;
                x.opts.bLoadData = false;
                x.opts.ProcessModules = [0 0 0];
                x.opts.ImportModules = [0 0 0 0];
            end
        end
    end
    
    % Try to find "dir" fields that were not found above
    if ~isfield(x.dir,'sourceStructure')
        x.dir.sourceStructure = '';
    end
    if ~isfield(x.dir,'studyPar')
        x.dir.studyPar = '';
    end
    if ~isfield(x.dir,'dataset_description')
        x.dir.dataset_description = '';
    end
    if ~isfield(x.dir,'dataPar')
        x.dir.dataPar = '';
    end

    % Recheck the JSON files (do they exist and which ones do)
    
    % Fallbacks
    bSourceStructure = false;
    bDatasetDescription = false;
    bDataPar = false;
    
    % Check if files exist
    if isfield(x,'dir')
        if isfield(x.dir,'sourceStructure')
            bSourceStructure = exist(x.dir.sourceStructure,'file');
        end
        if isfield(x.dir,'dataset_description')
            bDatasetDescription = exist(x.dir.dataset_description,'file');
        end
        if isfield(x.dir,'dataPar')
            bDataPar = exist(x.dir.dataPar,'file');
        end
    end
    
    if ~bSourceStructure && ~bDatasetDescription && ~bDataPar
        if x.opts.bImportData || x.opts.bProcessData
            fprintf('Neither the sourceStructure.json, dataset_description.json nor dataPar.json exist, ExploreASL will only be initialized...\n');
            % Check for really wrong input
            if ~isempty(regexp(x.opts.DatasetRoot, 'derivatives', 'once')) || ~isempty(regexp(x.opts.DatasetRoot, 'ExploreASL', 'once'))
                warning('Please do not provide the derivatives or ExploreASL folder. Use the study root directory instead...');
            end
        end
        x.opts.bProcessData = 0;
        x.opts.ProcessModules = [0 0 0];
    else
        % At least one of the JSON files exists
        
        % Dataset directory
        if strcmp(x.opts.dataParType,'directory')
            % Check if processing is turned off & there is a derivatives directory to be loaded
            if x.opts.bProcessData==0 && xASL_exist(fullfile(x.dir.DatasetRoot,'derivatives'))
                x.opts.bLoadData = true;
            end
        end
        
        % dataPar.json
        if strcmp(x.opts.dataParType,'dataParFile')
            % It is a dataPar.json, so do not run the BIDS import workflow
            if x.opts.bProcessData==0
                x.opts.bProcessData = 0; % Initialize & load but do not process
                x.opts.bLoadData = true;
            end
        end
        
        % sourceStructure.json
        if strcmp(x.opts.dataParType,'sourceStructure') || strcmp(x.opts.dataParType,'dataset_description')
            % It is a sourceStructure.json or dataset_description.json, so we run the import workflow
            if x.opts.bProcessData==0
                x.opts.bProcessData = 0; % Initialize & load but do not process
                x.opts.bLoadData = true;
            end
        end
        
    end

    % Check output
    if x.opts.bProcessData>0 && nargout==0
        warning('Data loading requested but no output structure defined');
        fprintf('%s\n', 'Try adding "x = " to the command to load data into the x structure');
    end
    
    % Try to catch unexpected inputs
    if strcmp(x.opts.dataParType,'unknown') && x.opts.bProcessData>0 && x.opts.bImportData==0
        fprintf('You are trying to process a dataset, without providing a dataPar.json file or running the import workflow...\n');
        x.opts.bProcessData = 0;
    end

end



