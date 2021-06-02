function [x] = xASL_init_DetermineRequiredPaths(x)
%xASL_init_DetermineRequiredPaths Check the BIDS dataset root for the metadata JSON files
%
% FORMAT: 
%   [x] = xASL_init_DetermineRequiredPaths(x)
%
% INPUT:
%   x             - ExploreASL x structure (REQUIRED, STRUCT)
%
% OUTPUT:
%   x             - ExploreASL x structure
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Check the BIDS dataset root for the metadata JSON files.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     n/a
%
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% Check the BIDS dataset root for the metadata JSON files

    % BIDS DatasetRoot directory
    x.dir.DatasetRoot = x.opts.DatasetRoot;
    x.opts.dataParType = 'directory';
    % Search for descriptive JSON files
    fileListSourceStructure = xASL_adm_GetFileList(x.dir.DatasetRoot, 'sourceStructure.*.json');
    fileListStudyPar = xASL_adm_GetFileList(x.dir.DatasetRoot, 'studyPar.*.json');
    fileListDataDescription = xASL_adm_GetFileList(fullfile(x.dir.DatasetRoot, 'rawdata'), 'dataset_description.json');
    % First try the derivatives folder
    fileListDataPar = xASL_adm_GetFileList(fullfile(x.dir.DatasetRoot, 'derivatives', 'ExploreASL'), 'dataPar*.json');
    % Check for x.D.ROOT
    if exist(fullfile(x.dir.DatasetRoot, 'derivatives', 'ExploreASL'),'dir')
        x.D.ROOT = fullfile(x.dir.DatasetRoot, 'derivatives', 'ExploreASL');
    end
    if isempty(fileListDataPar)
        % Derivatives maybe does not exist already, we'll try study root
        fileListDataPar = xASL_adm_GetFileList(x.dir.DatasetRoot, 'dataPar.*.json');
    end
    % Assign fields
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
        x.dir.dataPar = fileListDataPar{1};
    end

end


