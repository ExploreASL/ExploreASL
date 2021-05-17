function [x] = ExploreASL_Initialize(varargin)
%ExploreASL_Initialize Initializes ExploreASL
%
% FORMAT: 
%   [x] = ExploreASL_Initialize([DataParPath, ImportModules, ProcessModules, bPause, iWorker, nWorkers])
%
% INPUT:
%   varargin    - This script accepts the same arguments as ExploreASL_Master. Check out the definitions there.
%
% OUTPUT:
%   x           - struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: 
%
% This initialization wrapper initializes ExploreASL: managing paths, deployment, etc.
% Using the following initialization functions functions:
%
% xASL_init_DefinePaths            - manages folders for ExploreASL, and sets and creates study/data-folders
% xASL_init_Toolboxes              - initialization third-party toolboxes, e.g. SPM, dip_image (soon to be removed)
% xASL_init_VisualizationSettings  - defines visualization settings for
%                                    visual QC figure printing (type help xASL_init_VisualizationSettings for more information)
% xASL_init_DefineSets             - Define study subjects/parameters for this pipeline run
% xASL_init_PrintCheckSettings     - prints summarized data parameters and warnings
% xASL_init_FileSystem             - dirty initialization of common filenames used throughout the pipeline
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:
%
% Calling externally:             [x] = ExploreASL_Initialize('/MyDisk/MyStudy/DataPar.json');
% Debugging/initialization only:  [x] = ExploreASL_Initialize;
%
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% Admin

    % Define input parser
    p = inputParsing(varargin{:});

    % Convert parsed input
    parameters = ExploreASL_Initialize_convertParsedInput(p.Results);

    % Store parsed input
    x = ExploreASL_Initialize_storeParsedInput(parameters);
    
    % Initialize S substruct
    x.S = struct;
    
    % Check if the ExploreASL pipeline should be run or not
    if sum(x.ImportModules)>0
        x.bImportData = 1; % Importing data
        x.bReinitialize = true;
    else
        x.bImportData = 0; % Importing data
        x.bReinitialize = false;
    end
    if sum(x.ProcessModules)>0
        x.bProcessData = 1; % Loading & processing dataset
    else
        x.bProcessData = 0; % Only initialize ExploreASL functionality
    end
    if ~x.bProcessData && ~x.bImportData
        x.bReinitialize = false;
    end

    % Check if the DataParPath is a file or a directory
    SelectParFile = false; % Fallback
    if x.bProcessData
        % Checkout the "Proceed with Initialization" section
        if (isempty(x.DataParPath) || ~exist(x.DataParPath,'file'))
            SelectParFile = true; % If the DataParPath is either empty OR the file does not exist, we have to select it later on (if processing is turned on)
        end
    end
    

    %% -----------------------------------------------------------------------------
    %% Get ExploreASL path

    % Check if the current directory is the ExploreASL directory
    CurrCD = pwd;
    if exist(fullfile(CurrCD, 'ExploreASL_Master.m'), 'file')
        x.MyPath = CurrCD;
    end

    % Check whether MyPath is correct, otherwise obtain correct folder
    if ~isfield(x, 'MyPath')
        x.MyPath = '/DummyPath';
    end

    % Get the master script path
    MasterScriptPath = fullfile(x.MyPath, 'ExploreASL_Master.m');

    % Select the ExploreASL folder manually, if the script is not run in deployed mode
    if ~isdeployed
        if ~exist(MasterScriptPath,'file')
            pathstr = input('Provide foldername where ExploreASL is installed (format: \''PathExploreASL\''): ');
            if sum(pathstr==0) || ~exist(fullfile(pathstr,'ExploreASL_Master.m'),'file'), return; end
            x.MyPath = pathstr;
        end
    else
        % In deployed mode set the ExploreASL directory in the ctf archive
        [files,~] = spm_select('FPListRec',ctfroot,'ExploreASL_Master*'); % Find the path of the master files within the ctf archive
        curPathCTF = fileparts(files(1,:)); % Get the path
        x.MyPath = fullfile(curPathCTF); % curPathCTF = ExploreASL path

        BreakString = '==============================================================================================\n';
        fprintf(BreakString);
        fprintf('ctfroot:  %s\n', ctfroot);
        fprintf('x.MyPath: %s\n', x.MyPath);
        fprintf(BreakString);

    end

    % Go to ExploreASL folder
    cd(x.MyPath);


    %% Add ExploreASL paths
    if ~isdeployed
        addExploreASLDirectory(x.MyPath)
    end
    
    
    %% Check DataParPath
    [x] = ExploreASL_Initialize_checkDataParPath(x);
    
    
    % Give some feedback
    reportProcess = '';
    if x.bProcessData==0,        reportProcess = 'run the initialization';
    elseif x.bProcessData==1,    reportProcess = 'run the processing pipeline';
    elseif x.bProcessData==2,    reportProcess = 'load the dataset';
    end
    if x.bImportData==1,         reportImport = 'will run the import workflow and ';
    else,                       reportImport = '';
    end
    % Print feedback
    fprintf('ExploreASL %swill %s...\n',reportImport,reportProcess);
    
    %% Proceed with Initialization

    % Go to ExploreASL folder
    cd(x.MyPath);

    % Check if DataParFile needs to be loaded
    if x.bProcessData>0
        x = xASL_init_LoadDataParameterFile(x, x.DataParPath, SelectParFile);
    end


    %% Initialize general parameters
    x = xASL_init_DefineIndependentSettings(x); % these settings are data-independent

    x = xASL_init_DefineDataDependentSettings(x); % these settings depend on the data (e.g. which template to use)


    %% --------------------------------------------------------------------------------------------------------------------
    % Print logo
    BreakString = '==============================================================================================\n';

    LogoString = [...
    ' ________                      __                                 ______    ______   __        \n'...
    '/        |                    /  |                               /      \\  /      \\ /  |      \n'...
    '########/  __    __   ______  ## |  ______    ______    ______  /######  |/######  |## |      \n'...
    '## |__    /  \\  /  | /      \\ ## | /      \\  /      \\  /      \\ ## |__## |## \\__##/ ## |      \n'...
    '##    |   ##  \\/##/ /######  |## |/######  |/######  |/######  |##    ## |##      \\ ## |      \n'...
    '#####/     ##  ##<  ## |  ## |## |## |  ## |## |  ##/ ##    ## |######## | ######  |## |      \n'...
    '## |_____  /####  \\ ## |__## |## |## \\__## |## |      ########/ ## |  ## |/  \\__## |## |_____ \n'...
    '##       |/##/ ##  |##    ##/ ## |##    ##/ ## |      ##       |## |  ## |##    ##/ ##       |\n'...
    '########/ ##/   ##/ #######/  ##/  ######/  ##/        #######/ ##/   ##/  ######/  ########/ \n'...
    '                    ## |                                                                      \n'...
    '                    ## |                                                                      \n'...
    '                    ##/  \n'];

    fprintf([BreakString LogoString '\n']);

    %% Print chosen settings
    xASL_init_printSettings(x);


    %% Check permissions
    %xASL_adm_CheckPermissions(x.MyPath, false);


    %% Data-specific initialization
    fprintf('ExploreASL v%s initialized ... \n', x.Version);
    
    if x.bProcessData>0 && ~x.bImportData % Skip this step if we still need to run the import (first initialization)
        % Check if a root directory was defined
        if isempty(x.D.ROOT)
            error('No root/analysis/study folder defined');
        end

        % Fix a relative path
        if strcmp(x.D.ROOT(1), '.')
            cd(x.D.ROOT);
            x.D.ROOT = pwd;
        end

        % Define study subjects/parameters for this pipeline run
        x = xASL_init_DefineStudyData(x);


        % Remove lock dirs from previous runs, if ExploreASL is not running in parallel
        if x.nWorkers==1
            x = xASL_init_RemoveLockDirs(x);
        end

        % Define & print settings
        x = xASL_init_PrintCheckSettings(x);
        x = xASL_init_FileSystem(x);
    end

end


%% ==================================================================================
%% ==================================================================================
function [x] = xASL_init_RemoveLockDirs(x)
%xASL_init_RemoveLockDirs Remove 'lock-dir' if present from aborted previous run, for current subjects only


    % LockDir within 2 directories (e.g. T1, FLAIR or ASL)
    LockDir = fullfile(x.D.ROOT, 'lock');

    if exist(LockDir, 'dir')
        % fprintf('%s\n','Searching for locked previous ExploreASL image processing');
        LockDirFound = 0;
        LockDir = xASL_adm_FindByRegExp(fullfile(x.D.ROOT, 'lock'), {'(ASL|Structural|LongReg_T1)', x.subject_regexp, '.*module.*','^(locked)$'}, 'Match', 'Directories');
        if ~isempty(LockDir)
            warning('Locked folders were found, consider removing them before proceeding');
%                 for iL=1:length(LockDir)
%                     if isdir(LockDir{1})
%                         fprintf('%s\n',[LockDir{1} ' detected, removing']);
%                         rmdir(LockDir{1},'s');
%                     end
%                 end
%                 LockDirFound = 1;
        end

        % LockDir within 2 directories (e.g. DARTEL)
        LockDir = xASL_adm_FindByRegExp(fullfile(x.D.ROOT, 'lock'), {'(Population|DARTEL_T1)', '.*module.*','^(locked)$'}, 'Match','Directories');
        if ~isempty(LockDir)
            warning('Locked folders were found, consider removing them before proceeding');
%                 for iL=1:length(LockDir)
%                     fprintf('%s\n',[LockDir{1} ' detected, removing']);
%                     rmdir(LockDir{1},'s');
%                 end
%                 LockDirFound = 1;
        end

        if LockDirFound==0
            % fprintf('%s\n', 'No locked folders found from previous ExploreASL image processing');
        end
    end

end


%% -----------------------------------------------------------------------
%% Add ExploreASL Directory
function addExploreASLDirectory(MyPath)

    % Define paths (should be equal when loading data or only initializing)
    addpath(MyPath); % ExploreASL
    subfoldersToAdd = {...
        'Development', 'External', 'Functions', 'mex', 'Modules', 'Testing',...
        fullfile('Modules', 'SubModule_Import'), ...
        fullfile('Modules', 'SubModule_Structural'), ...
        fullfile('Modules', 'SubModule_ASL'), ...
        fullfile('Modules', 'SubModule_Population'), ...
        fullfile('External','isnear'), ...
        fullfile('External','DCMTK'), ...
        fullfile('Testing', 'UnitTests'), ...
        fullfile('External','ExploreQC'), ...
        fullfile('External','SPMmodified'), ...
        fullfile('External','SPMmodified','matlabbatch'),...
        fullfile('External','SPMmodified','xASL'),...
        fullfile('External','SPMmodified','toolbox','cat12'), ...
        fullfile('External','SPMmodified','toolbox','LST'), ...
        fullfile('External','SPMmodified','toolbox','OldNorm') ...
        genpath(fullfile('External','bids-matlab'))};
    % Iterate over subfolders which should be added
    for ii=1:length(subfoldersToAdd)
        addpath(fullfile(MyPath,subfoldersToAdd{ii}));
    end
    
end

%% -----------------------------------------------------------------------
%% Define input parser
function p = inputParsing(varargin)

    % Initialize input parser
    p = inputParser;
    
    % Define valid input variables
    validDataParPath = @(variable) ischar(variable) || isempty(variable);
    validImportModules = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable) || islogical(variable);
    validProcessModules = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable) || islogical(variable);
    validbPause = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable) || islogical(variable);
    validiWorker = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable);
    validnWorkers = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable);
    
    % Define defaults
    defaultDataParPath = [];
    defaultImportModules = [0 0 0 0];
    defaultProcessModules = [0 0 0];
    defaultbPause = 0;
    defaultiWorker = 1;
    defaultnWorkers = 1;
    
    % Add definitions to the input parser
    addOptional(p, 'DataParPath', defaultDataParPath, validDataParPath);
    addOptional(p, 'ImportModules', defaultImportModules, validImportModules);
    addOptional(p, 'ProcessModules', defaultProcessModules, validProcessModules);
    addOptional(p, 'bPause', defaultbPause, validbPause);
    addOptional(p, 'iWorker', defaultiWorker, validiWorker);
    addOptional(p, 'nWorkers', defaultnWorkers, validnWorkers);
    
    % Parse input
    parse(p,varargin{:});

end

%% -----------------------------------------------------------------------
%% Convert parsed input
function parameters = ExploreASL_Initialize_convertParsedInput(parameters)

    % Check if inputs are empty or chars
    if isempty(parameters.DataParPath),     parameters.DataParPath = '';                                    end
    if ischar(parameters.ImportModules),    parameters.ImportModules = str2num(parameters.ImportModules);   end
    if ischar(parameters.ProcessModules),   parameters.ProcessModules = str2num(parameters.ProcessModules); end
    if ischar(parameters.bPause),           parameters.bPause = str2num(parameters.bPause);                 end
    if ischar(parameters.iWorker),          parameters.iWorker = str2num(parameters.iWorker);               end
    if ischar(parameters.nWorkers),         parameters.nWorkers = str2num(parameters.nWorkers);             end
    
    % Check length of arrays (single digit input)
    if length(parameters.ImportModules)<4
        parameters.ImportModules = [parameters.ImportModules(1),...
                                    parameters.ImportModules(1),...
                                    parameters.ImportModules(1),...
                                    parameters.ImportModules(1)];
    end
    if length(parameters.ProcessModules)<3
        parameters.ProcessModules = [parameters.ProcessModules(1),...
                                     parameters.ProcessModules(1),...
                                     parameters.ProcessModules(1)];
    end
    
    % Make it impossible to set bPause to true in deployed mode
    if isdeployed
        parameters.bPause = 0;
    end

    if sum(parameters.ImportModules)~=0 || sum(parameters.ProcessModules)~=0
        % If import or processing is requested,
        % Do not allow continuing without a valid JSON path as input
        [~, ~, Fext] = fileparts(parameters.DataParPath);
        if isempty(Fext) || ~strcmpi(Fext, '.json')
            error('Invalid path, first argument should be a path to a JSON-file');
        end
    end
    
end

%% -----------------------------------------------------------------------
%% Store parsed input
function x = ExploreASL_Initialize_storeParsedInput(parameters)

    % Store input
    x.DataParPath = parameters.DataParPath;
    x.ImportModules = parameters.ImportModules;
    x.ProcessModules = parameters.ProcessModules;
    x.bPause = parameters.bPause;
    x.iWorker = parameters.iWorker;
    x.nWorkers = parameters.nWorkers;
    
end


%% -----------------------------------------------------------------------
function [x] = ExploreASL_Initialize_checkDataParPath(x)


	% Check if the DataParPath is a directory (NEW - ASL BIDS)
	if exist(x.DataParPath,'dir'))
		% ASL-BIDS studyRoot directory
		x.StudyRoot = x.DataParPath;
	elseif exist(x.DataParPath,'file'))
		% Input is either a sourceStructure.json, dataset_description.json or dataPar.json
		fileListSourceStructure = xASL_adm_GetFileList(x.StudyRoot, 'sourceStructure*.json');
		fileListStudyPar = xASL_adm_GetFileList(x.StudyRoot, 'studyPar*.json');
		fileListDataPar = xASL_adm_GetFileList(x.StudyRoot, 'dataPar*.json');
		% ...
	end


	% Check if pipeline should be run, but there is no DataParPath
    if x.bProcessData && (~isfield(x,'DataParPath') || ~exist(x.DataParPath,'file'))
        x.DataParPath = input('Please insert the path to your DataParFile: ');
    end
    
    % Read file
    if exist(x.DataParPath,'file')==2
        SelectParFile = false; % Does not need to be inserted a second time
        jsonContent = spm_jsonread(x.DataParPath);
        if isfield(jsonContent,'x')
            x.dataParType = 'dataParFile';
        else
            x.dataParType = 'sourceStructure';
        end
    else
        x.dataParType = 'unknown';
    end

    % Recheck the DataPar/sourceStructure file, which is possibly not a file or does not exist
    if ~exist(x.DataParPath,'file')
        if x.bImportData || x.bProcessData
            fprintf('DataPar file does not exist, ExploreASL will only be initialized...\n');
        end
        x.bProcessData = 0;
        x.ProcessModules = [0 0 0];
    else % DataPar file/folder exists
        if strcmp(x.dataParType,'dataParFile') % It is a dataParFile, so do not run the BIDS import workflow
            if x.bProcessData==0 || x.bProcessData==2
                x.bProcessData = 2; % Initialize & load but do not process
                x.bReinitialize = false; % Do not reinitialize if we only load the data
            end
        end
    end

    % Check output
    if x.bProcessData>0 && nargout==0
        warning('Data loading requested but no output structure defined');
        fprintf('%s\n', 'Try adding "x = " to the command to load data into the x structure');
    end
    
    % Try to catch unexpected inputs
    if strcmp(x.dataParType,'unknown') && x.bProcessData>0 && x.bImportData==0
        fprintf('You are trying to process a dataset, without providing a DataPar file or running the import workflow...\n');
        x.bProcessData = 0;
    end


end


