function [x] = xASL_init_LoadDataPar(x)
%xASL_init_LoadDataPar Load data parameter file
%
% FORMAT: [x] = xASL_init_LoadDataParameterFile(x)
%
% INPUT:
%   x       - ExploreASL x structure (STRUCT, REQUIRED)
%
% OUTPUT:
%   x       - ExploreASL x structure (STRUCT)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Load data parameter file (which contains settings for running the pipeline, 
%              so should actually be called something like xASLProcessingSettings.json
%
% EXAMPLE:     This is part of the initialization workflow. Check out the usage there.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:  n/a
% __________________________________
% Copyright (c) 2015-2024 ExploreASL
    
    

    %% Generate warning for incompatibility with old dataPar.m
    [~, ~, Dext] = fileparts(x.dir.dataPar);
    if strcmp(Dext,'.m')
        warning('No .m file backwards compatibility starting v1.10.0...');
    elseif strcmp(Dext,'.mat')
        warning('No .mat file backwards compatibility starting v1.10.0...');
    end


    %% Choose the dataPar location

    bUseRoot = false;
    bUseDerivatives = false;

    % Check dataPar inside the root folder /
    listDataparRoot = xASL_adm_GetFileList(x.dir.DatasetRoot, '(?i)(^dataPar.*\.json$)', 'FPList', [], 0);
    if ~isempty(listDataparRoot)
        bUseRoot = true;
    end   

    % Check dataPar inside the rawdata folder /rawdata (this one we ignore, ExploreASL settings shouldn't belong to a raw BIDS dataset)
    listDataparRawdata = xASL_adm_GetFileList(x.dir.RawData, '(?i)(^dataPar.*\.json$)', 'FPList', [], 0);
    if ~isempty(listDataparRawdata)
        warning([xASL_num2str(length(listDataparRawdata)) ' dataPar.json (or similar) file(s) found in ' x.dir.RawData ', will try this but this is not the appropriate location']);
        bUseRawdata = true;
    end

    % Check dataPar inside the processing folder /derivatives/ExploreASL
    fListDataParLegacy = xASL_adm_GetFileList(x.dir.xASLDerivatives, '(?i)(^dataPar.*\.json$)', 'FPList', [], 0);
    if ~isempty(fListDataParLegacy)
        bUseDerivatives = true;
    end

    % Choose folder & filelist
    if bUseDerivatives && bUseRoot
        fprintf('%s\n', 'dataPar.json (or similar) file(s) both in /derivatives/ExploreASL and / root folder, ignoring the last');
        bUseRoot = false;
    end
    if bUseDerivatives
        fFolder = x.dir.xASLDerivatives;
        listDatapar = fListDataParLegacy;
        fprintf('%s\n', ['dataPar.json found in ' fFolder]);
    elseif bUseRoot
        fFolder = x.dir.DatasetRoot;
        listDatapar = listDataparRoot;
        fprintf('%s\n', ['dataPar.json found in ' fFolder]);
    elseif bUseRawdata
        fFolder = x.dir.RawData;
        listDatapar = listDataparRawdata;
        fprintf('%s\n', ['dataPar.json found in ' fFolder]);
    else
        listDatapar = {};
    end


    %% Load pre-existing dataPar
    if length(listDatapar)>1
        fprintf('Warning: multiple dataPar.json files found, using the first\n');
    elseif isempty(listDatapar)
        % Create default if no dataPar was provided
        fprintf('No dataPar.json provided, will create default version\n');
        dataPar = struct();
    else
        % dataPar was provided
        dataPar.x = xASL_io_ReadDataPar(listDatapar{1});
    end


    %% Populate dataPar with missing parameters

    % Fills in important information in the dataPar if missing
    if ~isfield(dataPar, 'x')
        % Add x field
        dataPar.x = struct;
    end
    % Dataset fields
    if ~isfield(dataPar.x,'dataset')
        dataPar.x.dataset = struct;
    end
    % Check for settings fields
    if ~isfield(dataPar.x,'settings')
        dataPar.x.settings = struct;
    end
    % Check for quality field
    if ~isfield(dataPar.x.settings,'Quality')
        dataPar.x.settings.Quality = true;
    end
    % Check for DELETETEMP field
    if ~isfield(dataPar.x.settings,'DELETETEMP')
        dataPar.x.settings.DELETETEMP = true;
    end

    %% Write final dataPar
    x.dir.dataPar = fullfile(x.dir.xASLDerivatives, 'dataPar.json');
    xASL_io_WriteJson(x.dir.dataPar, dataPar);

    
end