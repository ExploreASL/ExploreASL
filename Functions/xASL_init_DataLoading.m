function [x] = xASL_init_DataLoading(x)
%xASL_init_DataLoading Load dataset by adding relevant fields to xASL x struct
%
% FORMAT: [x] = xASL_init_DataLoading(x)
% 
% INPUT:
%   x          - ExploreASL x structure (STRUCT, REQUIRED)
%
% OUTPUT:
%   x          - ExploreASL x structure
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Load dataset by adding relevant fields to xASL x struct.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        [x] = xASL_init_DataLoading(x);
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% Data loading
    if ~isfield(x,'dataset')
        x.dataset = struct;
    end
    
    % Make sure that the dataPar.json definitely exists if we load the dataset
    if x.opts.bLoadData
        if isfield(x,'dir') && isfield(x.dir,'dataPar')
            if isempty(x.dir.dataPar)
                if ~x.opts.bImportData
                    warning('You are trying to load a dataset but there is no dataPar JSON file...');
                end
                x.opts.bLoadData = false;
            end
        else
            x.opts.bLoadData = false;
        end
    end
    
    % Go to ExploreASL folder
    cd(x.opts.MyPath);

    % Check if DataParFile needs to be loaded
    if x.opts.bProcessData || x.opts.bLoadData
        if ~isempty(x.dir.dataPar)
            x = xASL_init_LoadDataParameterFile(x, x.dir.dataPar, x.settings.SelectParFile);
        else
            fprintf('No dataPar.json provided...\n');
            if x.opts.bLoadData
                fprintf('Dataset can not be loaded...\n');
                x.opts.bLoadData = false;
            end
        end
    end
    
    % These settings depend on the data (e.g. which template to use)
    x = xASL_init_DefineDataDependentSettings(x);
    
    % Check if a "loadable" dataset exists (xASL_exist(x.dir.dataPar) is only there for backwards compatibility)
    loadableDataset = (isfield(x.dir,'DatasetRoot') && xASL_exist(fullfile(x.dir.DatasetRoot,'derivatives'),'dir')) ...
                   || (~isempty(x.dir.dataPar) && xASL_exist(x.dir.dataPar,'file'));
    
    % Check if data loading should be executed first
    if x.opts.bLoadData && loadableDataset
        % Check if a root directory was defined
        if ~isfield(x.D,'ROOT') || isempty(x.D.ROOT)
            error('No root folder defined');
        end

        % Fix a relative path
        if strcmp(x.D.ROOT(1), '.')
            cd(x.D.ROOT);
            x.D.ROOT = pwd;
        end

        % Define study subjects/parameters for this pipeline run
        x = xASL_init_DefineStudyData(x);

        % Remove lock dirs from previous runs, if ExploreASL is not running in parallel
        if x.opts.nWorkers==1
            x = xASL_init_RemoveLockDirs(x);
        end

        % Define & print settings
        x = xASL_init_PrintCheckSettings(x);
    elseif x.opts.bLoadData && ~loadableDataset
        % This warning is also printed if a user tries to "only load" a dataset with a descriptive JSON file. 
        % Since this behavior will be discontinued (only directories from now on), I do not see a problem with this for now.
        warning('Dataset can not be loaded, there is no derivatives directory, try to run the import first...');
    end


end


