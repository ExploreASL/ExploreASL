function [x] = ExploreASL_Initialize(varargin)
%ExploreASL_Initialize Initializes ExploreASL
%
% FORMAT: 
%   [x] = ExploreASL_Initialize([StudyRoot, ImportModules, ProcessModules, bPause, iWorker, nWorkers])
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
    x = ExploreASL_Initialize_GetBooleansImportProcess(x);

    % Check if the StudyRoot is a file or a directory
    SelectParFile = false; % Fallback
    if x.opts.bProcessData
        % Checkout the "Proceed with Initialization" section
        if (isempty(x.opts.StudyRoot) || (~exist(x.opts.StudyRoot,'file') && ~exist(x.opts.StudyRoot,'dir')))
            SelectParFile = true; % If the StudyRoot is either empty OR the file does not exist, we have to select it later on (if processing is turned on)
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
        % Check if we can get the path from the ExploreASL_Initialize path
        initializePath = fileparts(mfilename('fullpath'));
        if ~isempty(regexp(initializePath,'ExploreASL$', 'once'))
            x.MyPath = initializePath;
        else
            x.MyPath = '/DummyPath';
        end
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
    
    
    %% Check StudyRoot
    [x, SelectParFile] = ExploreASL_Initialize_checkStudyRoot(x, SelectParFile);
    
    % Give some feedback
    ExploreASL_Initialize_basicFeedback(x);
    
    %% Proceed with Initialization

    % Go to ExploreASL folder
    cd(x.MyPath);

    % Check if DataParFile needs to be loaded
    if x.opts.bProcessData>0
        if ~isempty(x.dir.dataPar)
            x = xASL_init_LoadDataParameterFile(x, x.dir.dataPar, SelectParFile);
        else
            fprintf('No dataPar.json provided...\n');
        end
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
    
    if x.opts.bProcessData>0 && ~x.opts.bImportData % Skip this step if we still need to run the import (first initialization)
        % Check if a root directory was defined
        if ~isfield(x.D,'ROOT') || isempty(x.D.ROOT)
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
        if x.opts.nWorkers==1
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
        end

        % LockDir within 2 directories (e.g. DARTEL)
        LockDir = xASL_adm_FindByRegExp(fullfile(x.D.ROOT, 'lock'), {'(Population|DARTEL_T1)', '.*module.*','^(locked)$'}, 'Match','Directories');
        if ~isempty(LockDir)
            warning('Locked folders were found, consider removing them before proceeding');
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
    validStudyRoot = @(variable) ischar(variable) || isempty(variable);
    validImportModules = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable) || islogical(variable);
    validProcessModules = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable) || islogical(variable);
    validbPause = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable) || islogical(variable);
    validiWorker = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable);
    validnWorkers = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable);
    
    % Define defaults
    defaultStudyRoot = [];
    defaultImportModules = [0 0 0 0];
    defaultProcessModules = [0 0 0];
    defaultbPause = 0;
    defaultiWorker = 1;
    defaultnWorkers = 1;
    
    % Add definitions to the input parser
    addOptional(p, 'StudyRoot', defaultStudyRoot, validStudyRoot);
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
    if isempty(parameters.StudyRoot),     parameters.StudyRoot = '';                                    end
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
    
    % Check nWorkers
    if numel(parameters.nWorkers)>1
        warning('nWorkers is supposed to be an integer number, resetting to 1...');
        parameters.nWorkers = 1;
    end
    
end

%% -----------------------------------------------------------------------
%% Store parsed input
function x = ExploreASL_Initialize_storeParsedInput(parameters)

    % Store input options
    x.opts.StudyRoot = parameters.StudyRoot;
    x.opts.ImportModules = parameters.ImportModules;
    x.opts.ProcessModules = parameters.ProcessModules;
    x.opts.bPause = parameters.bPause;
    x.opts.iWorker = parameters.iWorker;
    x.opts.nWorkers = parameters.nWorkers;
    
end


%% -----------------------------------------------------------------------
function [x, SelectParFile] = ExploreASL_Initialize_checkStudyRoot(x, SelectParFile)

    % Check if the StudyRoot is a directory (NEW - ASL BIDS)
    x.dataParType = 'unknown'; % Fallback
    % Create directory field if it doesn't exist already
    if ~isfield(x, 'dir')
        x.dir = struct;
    end
    if x.opts.bImportData || x.opts.bProcessData
        if exist(x.opts.StudyRoot,'dir')
            % ASL-BIDS studyRoot directory
            x.dir.StudyRoot = x.opts.StudyRoot;
            x.dataParType = 'directory';
            % Search for descriptive JSON files
            fileListSourceStructure = xASL_adm_GetFileList(x.dir.StudyRoot, 'sourceStructure.*.json');
            fileListStudyPar = xASL_adm_GetFileList(x.dir.StudyRoot, 'studyPar.*.json');
            fileListDataDescription = xASL_adm_GetFileList(fullfile(x.dir.StudyRoot, 'rawdata'), 'dataset_description.json');
            % First try the derivatives folder
            fileListDataPar = xASL_adm_GetFileList(fullfile(x.dir.StudyRoot, 'derivatives', 'ExploreASL'), 'dataPar*.json');
            if isempty(fileListDataPar)
                % Derivatives maybe does not exist already, we'll try study root
                fileListDataPar = xASL_adm_GetFileList(x.dir.StudyRoot, 'dataPar.*.json');
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
        elseif exist(x.opts.StudyRoot,'file')
            % Temporary functionality, this will lead to an error starting v2.0.0
            [x, SelectParFile] = ExploreASL_Initialize_checkStudyRoot_invalid_starting_2_0(x);
        else
            if x.opts.bProcessData || x.opts.bImportData
                if ~isdeployed
                    x.opts.StudyRoot = input('Please insert the path to your study directory: ');
                else
                    error('Study directory does not exist...');
                end
                % Immediately check the input
                if ~exist(x.opts.StudyRoot, 'dir')
                    warning('This study directory does not exist, ExploreASL will only be initialized...');
                    x.opts.bProcessData = 0;
                    x.opts.bImportData = 0;
                    x.opts.bReinitialize = 0;
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
        end
        x.opts.bProcessData = 0;
        x.opts.ProcessModules = [0 0 0];
    else
        % At least one of the JSON files exists
        
        if strcmp(x.dataParType,'dataParFile')
            % It is a dataPar.json, so do not run the BIDS import workflow
            if x.opts.bProcessData==0 || x.opts.bProcessData==2
                x.opts.bProcessData = 2; % Initialize & load but do not process
                x.bReinitialize = false; % Do not reinitialize if we only load the data
            end
        end
        
        if strcmp(x.dataParType,'sourceStructure') || strcmp(x.dataParType,'dataset_description')
            % It is a sourceStructure.json or dataset_description.json, so we run the import workflow
            if x.opts.bProcessData==0 || x.opts.bProcessData==2
                x.opts.bProcessData = 2; % Initialize & load but do not process
                x.bReinitialize = true; % Do not reinitialize if we only load the data
            end
        end
        
    end

    % Check output
    if x.opts.bProcessData>0 && nargout==0
        warning('Data loading requested but no output structure defined');
        fprintf('%s\n', 'Try adding "x = " to the command to load data into the x structure');
    end
    
    % Try to catch unexpected inputs
    if strcmp(x.dataParType,'unknown') && x.opts.bProcessData>0 && x.opts.bImportData==0
        fprintf('You are trying to process a dataset, without providing a dataPar.json file or running the import workflow...\n');
        x.opts.bProcessData = 0;
    end


end


%% -----------------------------------------------------------------------
% Check if the ExploreASL pipeline should be run or not
function [x] = ExploreASL_Initialize_GetBooleansImportProcess(x)

    if sum(x.opts.ImportModules)>0
        x.opts.bImportData = 1; % Importing data
        x.opts.bReinitialize = true;
    else
        x.opts.bImportData = 0; % Importing data
        x.opts.bReinitialize = false;
    end
    if sum(x.opts.ProcessModules)>0
        x.opts.bProcessData = 1; % Loading & processing dataset
    else
        x.opts.bProcessData = 0; % Only initialize ExploreASL functionality
    end
    if ~x.opts.bProcessData && ~x.opts.bImportData
        x.opts.bReinitialize = false;
    end

end


%% -----------------------------------------------------------------------
% Give some feedback
function ExploreASL_Initialize_basicFeedback(x)

    % Report string
    reportProcess = '';
    if x.opts.bProcessData==0
        reportProcess = 'run the initialization';
    elseif x.opts.bProcessData==1
        reportProcess = 'run the processing pipeline';
    elseif x.opts.bProcessData==2
        reportProcess = 'load the dataset';
    end
    if x.opts.bImportData==1
        reportImport = 'will run the import workflow and ';
    else
        reportImport = '';
    end
    % Print feedback
    fprintf('ExploreASL %swill %s...\n',reportImport,reportProcess);

end


%% -----------------------------------------------------------------------
function [x, SelectParFile] = ExploreASL_Initialize_checkStudyRoot_invalid_starting_2_0(x)

    % Input is either a sourceStructure.json, dataset_description.json or dataPar.json
    warning('You provided a descriptive JSON file. We recommend to use the study root folder instead...');
    SelectParFile = false; % Does not need to be inserted a second time
    [~, ~, extensionJSON] = fileparts(x.opts.StudyRoot);
    if strcmp(extensionJSON,'.json') || strcmp(extensionJSON,'.JSON')
        % Try to find out type by name
        if ~isempty(regexp(x.opts.StudyRoot, 'sourceStructure', 'once'))
            x.dataParType = 'sourceStructure';
            x.dir.sourceStructure = x.opts.StudyRoot;
        elseif ~isempty(regexp(x.opts.StudyRoot, 'studyPar', 'once'))
            x.dataParType = 'studyPar';
            x.dir.studyPar = x.opts.StudyRoot;
            warning('You provided the studyPar.json, which should never be the input...');
        elseif ~isempty(regexp(x.opts.StudyRoot, 'dataset_description', 'once'))
            x.dataParType = 'dataset_description';
            x.dir.dataset_description = x.opts.StudyRoot;
        elseif ~isempty(regexp(x.opts.StudyRoot, 'dataPar', 'once'))
            x.dataParType = 'dataParFile';
            x.dir.dataPar = x.opts.StudyRoot;
        else
            % No files with correct names found
            error('No matching JSON files found...');
        end
    end
    
    % Try to find study root from existing files
    if isfield(x.dir,'sourceStructure')
        % We expect the sourceStructure.json to be within the study root folder
        [x.dir.StudyRoot, ~] = fileparts(x.dir.sourceStructure);
        x.opts.StudyRoot = x.dir.StudyRoot;
    end
    if isfield(x.dir,'studyPar')
        % We expect the sourceStructure.json to be within the study root folder
        [x.dir.StudyRoot, ~] = fileparts(x.dir.studyPar);
        x.opts.StudyRoot = x.dir.StudyRoot;
    end
    if isfield(x.dir,'dataset_description')
        % We expect the dataset_description.json to be within the rawdata folder
        [rawdataFolder, ~] = fileparts(x.dir.dataset_description);
        x.dir.StudyRoot = fileparts(rawdataFolder);
    end
    if isfield(x.dir,'dataParFile')
        % We expect the dataPar.json to be within the study root folder
        [x.dir.StudyRoot, ~] = fileparts(x.dir.dataParFile);
    end
    
    % Recheck for other files if studyRoot is known now
    if isfield(x, 'dir')
        if isfield(x.dir, 'StudyRoot')
            fileListSourceStructure = xASL_adm_GetFileList(x.dir.StudyRoot, 'sourceStructure*.json');
            fileListStudyPar = xASL_adm_GetFileList(x.dir.StudyRoot, 'studyPar*.json');
            fileListDataDescription = xASL_adm_GetFileList(fullfile(x.dir.StudyRoot, 'rawdata'), 'dataset_description.json');
            fileListDataPar = xASL_adm_GetFileList(fullfile(x.dir.StudyRoot, 'derivatives', 'ExploreASL'), 'dataPar*.json');
            if isempty(fileListDataPar)
                % Derivatives maybe does not exist already, we'll try study root
                fileListDataPar = xASL_adm_GetFileList(x.dir.StudyRoot, 'dataPar*.json');
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
                x.dir.dataPar = fileListDataPar{1};
            end
        end
    end

end


