function [x, SelectParFile] = xASL_init_checkDatasetRoot_invalid_starting_2_0(x)
%xASL_init_checkDatasetRoot_invalid_starting_2_0 Script for backwards compatibility. 
% This functionality allows loading JSON files instead of the BIDS dataset root.
%
% FORMAT: 
%   [x, SelectParFile] = xASL_init_checkDatasetRoot_invalid_starting_2_0(x)
%
% INPUT:
%   x             - ExploreASL x structure (REQUIRED, STRUCT)
%
% OUTPUT:
%   x             - ExploreASL x structure
%   SelectParFile - Variable to determine if we need to ask the user another time for the DatasetRoot
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Script for backwards compatibility. This functionality allows loading JSON files instead of the BIDS dataset root.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     n/a
%
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% Find out which descriptive JSON file was provided, then try set up ExploreASL using it

    % Input is either a sourceStructure.json, dataset_description.json or dataPar.json
    warning('You provided a descriptive JSON file. We recommend to use the dataset root folder instead...');
    SelectParFile = false; % Does not need to be inserted a second time
    [~, ~, extensionJSON] = fileparts(x.opts.DatasetRoot);
    if strcmp(extensionJSON,'.json') || strcmp(extensionJSON,'.JSON')
        % Try to find out type by name
        if ~isempty(regexp(x.opts.DatasetRoot, 'sourceStructure', 'once'))
            x.opts.dataParType = 'sourceStructure';
            x.dir.sourceStructure = x.opts.DatasetRoot;
        elseif ~isempty(regexp(x.opts.DatasetRoot, 'studyPar', 'once'))
            x.opts.dataParType = 'studyPar';
            x.dir.studyPar = x.opts.DatasetRoot;
            warning('You provided the studyPar.json, which should never be the input...');
        elseif ~isempty(regexp(x.opts.DatasetRoot, 'dataset_description', 'once'))
            x.opts.dataParType = 'dataset_description';
            x.dir.dataset_description = x.opts.DatasetRoot;
        elseif ~isempty(regexpi(x.opts.DatasetRoot, 'dataPar', 'once'))
            x.opts.dataParType = 'dataParFile';
            x.dir.dataPar = x.opts.DatasetRoot;
        else
            % No files with correct prefix in the filename found
            error('Could not identify the type of the provided JSON file:\n %s\n The name should start with sourceStructure, studyPar, dataset_description or dataPar',x.opts.DatasetRoot);
        end
    end
    
    % Try to find study root from existing files
    if isfield(x.dir,'sourceStructure')
        % We expect the sourceStructure.json to be within the study root folder
        [x.dir.DatasetRoot, ~] = fileparts(x.dir.sourceStructure);
        x.opts.DatasetRoot = x.dir.DatasetRoot;
    end
    if isfield(x.dir,'studyPar')
        % We expect the sourceStructure.json to be within the study root folder
        [x.dir.DatasetRoot, ~] = fileparts(x.dir.studyPar);
        x.opts.DatasetRoot = x.dir.DatasetRoot;
    end
    if isfield(x.dir,'dataset_description')
        % We expect the dataset_description.json to be within the rawdata folder
        [rawdataFolder, ~] = fileparts(x.dir.dataset_description);
        x.dir.DatasetRoot = fileparts(rawdataFolder);
    end
    if isfield(x.dir,'dataPar')
        % We expect the dataPar.json to be within the study root folder
        [x.dir.DatasetRoot, ~] = fileparts(x.dir.dataPar);
    end
    
    % Recheck for other files if DatasetRoot is known now
    if isfield(x, 'dir')
        if isfield(x.dir, 'DatasetRoot')
            fileListSourceStructure = xASL_adm_GetFileList(x.dir.DatasetRoot, 'sourceStructure*.json');
            fileListStudyPar = xASL_adm_GetFileList(x.dir.DatasetRoot, 'studyPar*.json');
            fileListDataDescription = xASL_adm_GetFileList(fullfile(x.dir.DatasetRoot, 'rawdata'), 'dataset_description.json');
            fileListDataPar = xASL_adm_GetFileList(fullfile(x.dir.DatasetRoot, 'derivatives', 'ExploreASL'), 'dataPar*.json');
            if isempty(fileListDataPar)
                % Derivatives maybe does not exist already, we'll try study root
                fileListDataPar = xASL_adm_GetFileList(x.dir.DatasetRoot, 'dataPar*.json');
            end
            if ~isempty(fileListSourceStructure)
                x.dir.sourceStructure = fileListSourceStructure{1};
            end
            if ~isempty(fileListStudyPar)
                x.dir.studyPar = fileListStudyPar{1};
            end
            if ~isempty(fileListDataDescription)
                x.dir.dataset_description = fileListDataDescription{1};
            end
            if ~isempty(fileListDataPar)
				if length(fileListDataPar) > 1
					warning('Multiple dataPar*.jsons exist. Using the first: %s\n',fileListDataPar{1});
				end
                x.dir.dataPar = fileListDataPar{1};
            end
        end
    end

end

