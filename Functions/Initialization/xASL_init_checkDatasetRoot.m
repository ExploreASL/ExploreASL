function [x] = xASL_init_checkDatasetRoot(x)
%xASL_init_checkDatasetRoot Check the ExploreASL parameter DatasetRoot
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
% Copyright (c) 2015-2022 ExploreASL

    %% Check the ExploreASL parameter "DatasetRoot"
    
    % Default
    x.opts.dataParType = 'unknown';
    
    % Create directory field if it doesn't exist already
    if ~isfield(x, 'dir')
        x.dir = struct;
    end

    % Check if user correctly inserted a dataset root directory
    if exist(x.opts.DatasetRoot,'dir')
        [x] = xASL_init_DetermineRequiredPaths(x);
    elseif exist(x.opts.DatasetRoot,'file')
        % Files are no longer supported for the dataset root directory
        error('You provided a descriptive JSON or another file. We only allow paths to the dataset root directory...');
    else
        % Check if the user inserted a directory or file which does not exist
        if ~isempty(x.opts.DatasetRoot)
            warning('Dataset root directory does not exist...');
        end
        if x.opts.bProcessData || x.opts.bImportData
            % Give back a warning that the user tried to import or process but neither a correct dataset root nor a dataPar.json that exists was used
            warning('You are trying to import or process a dataset, but the input parameters are not correct. ExploreASL will only be initialized...');
            x.opts.bProcessData = 0;
            x.opts.bImportData = 0;
            x.opts.bLoadData = false;
            x.opts.bProcess = [0 0 0];
            x.opts.bImport = [0 0 0];
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
            % Check for wrong input
            [~, DatasetDir] = fileparts(x.opts.DatasetRoot);
            if strcmp(DatasetDir, 'derivatives') || strcmp(DatasetDir, 'ExploreASL')
                warning('Please do not provide the derivatives or ExploreASL folder. Use the dataset root directory instead...');
            end
        end
        x.opts.bProcessData = 0;
        x.opts.bProcess = [0 0 0];
    else
        % At least one of the JSON files exists
        
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


