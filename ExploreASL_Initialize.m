function [x] = ExploreASL_Initialize(varargin)
%ExploreASL_Initialize Initializes ExploreASL
%
% FORMAT: 
%   [x] = ExploreASL_Initialize([DatasetRoot, ImportModules, ProcessModules, bPause, iWorker, nWorkers])
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
% This initialization workflow initializes ExploreASL. The overall workflow is shown below:
% 
% 1. Admin
% - Input parsing (inputParsing, ExploreASL_Initialize_convertParsedInput, ExploreASL_Initialize_storeParsedInput)
% - Initialize substructs of the ExploreASL x structure (ExploreASL_Initialize_SubStructs)
% - Check input parameters and determine pipeline related booleans (ExploreASL_Initialize_GetBooleansImportProcess)
%
% 2. Get ExploreASL path
% - Check if the current directory is the ExploreASL directory
% - Check whether MyPath is correct, otherwise obtain correct folder
% - If something went wrong, select the ExploreASL folder manually
%
% 3. Add ExploreASL paths
%
% 4. Check DatasetRoot
% - Check the provided DatasetRoot parameter (is it a path? does it exist? is it a directory?)
% - Print some basic feedback regarding the chosen ExploreASL settings
%
% 5. Check ExploreASL parameter file & general settings
% - Check if the dataPar.json exists and if it can be loaded or not
% - Initialize the data independent & dependent settings
%
% 6. Print logo & settings
%
% 7. Data-specific initialization
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:
%
% Calling externally:             [x] = ExploreASL_Initialize('/MyDisk/MyStudy');
% Debugging/initialization only:  [x] = ExploreASL_Initialize;
%
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% 1. Admin

    % Define input parser
    p = inputParsing(varargin{:});

    % Convert parsed input
    parameters = ExploreASL_Initialize_convertParsedInput(p.Results);

    % Store parsed input
    x = ExploreASL_Initialize_storeParsedInput(parameters);
    
    % Initialize substructs
    x = ExploreASL_Initialize_SubStructs(x);
    
    % Check if the ExploreASL pipeline should be run or not
    x = ExploreASL_Initialize_GetBooleansImportProcess(x);

    % Check if the DatasetRoot is a file or a directory
    SelectParFile = false; % Fallback
    if x.opts.bProcessData
        % Checkout the "Proceed with Initialization" section
        if (isempty(x.opts.DatasetRoot) || (~exist(x.opts.DatasetRoot,'file') && ~exist(x.opts.DatasetRoot,'dir')))
            SelectParFile = true; % If the DatasetRoot is either empty OR the file does not exist, we have to select it later on (if processing is turned on)
        end
    end
    
    
    %% 2. Get ExploreASL path

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


    %% 3. Add ExploreASL paths
    if ~isdeployed
        addExploreASLDirectory(x.MyPath)
    end
    
    
    %% 4. Check DatasetRoot
    [x, SelectParFile] = xASL_init_checkDatasetRoot(x, SelectParFile);
    
    % Give some feedback
    ExploreASL_Initialize_basicFeedback(x);
    
    %% 5. Check ExploreASL parameter file & general settings

    % Go to ExploreASL folder
    cd(x.MyPath);

    % Check if DataParFile needs to be loaded
    if x.opts.bProcessData || x.opts.bOnlyLoad
        if ~isempty(x.dir.dataPar)
            x = xASL_init_LoadDataParameterFile(x, x.dir.dataPar, SelectParFile);
        else
            fprintf('No dataPar.json provided...\n');
            if x.opts.bOnlyLoad
                fprintf('Dataset can not be loaded...\n');
                x.opts.bOnlyLoad = 0;
            end
        end
    end

    % Initialize general settings
    x = xASL_init_DefineIndependentSettings(x); % these settings are data-independent

    x = xASL_init_DefineDataDependentSettings(x); % these settings depend on the data (e.g. which template to use)


    %% 6. Print logo & settings
    
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

    % Print chosen settings
    xASL_init_printSettings(x);


    %% Check permissions
    %xASL_adm_CheckPermissions(x.MyPath, false);


    %% 7. Data-specific initialization
    fprintf('ExploreASL v%s initialized ... \n', x.Version);
    
    if (x.opts.bProcessData || x.opts.bOnlyLoad) && ~x.opts.bImportData % Skip this step if we still need to run the import (first initialization)
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
        x = xASL_init_FileSystem(x);
    end

end


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
    validDatasetRoot = @(variable) ischar(variable) || isempty(variable);
    validImportModules = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable) || islogical(variable);
    validProcessModules = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable) || islogical(variable);
    validbPause = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable) || islogical(variable);
    validiWorker = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable);
    validnWorkers = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable);
    
    % Define defaults
    defaultDatasetRoot = [];
    defaultImportModules = [0 0 0 0];
    defaultProcessModules = [0 0 0];
    defaultbPause = 0;
    defaultiWorker = 1;
    defaultnWorkers = 1;
    
    % Add definitions to the input parser
    addOptional(p, 'DatasetRoot', defaultDatasetRoot, validDatasetRoot);
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
    if isempty(parameters.DatasetRoot),     parameters.DatasetRoot = '';                                    end
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
    x.opts.DatasetRoot = parameters.DatasetRoot;
    x.opts.ImportModules = parameters.ImportModules;
    x.opts.ProcessModules = parameters.ProcessModules;
    x.opts.bPause = parameters.bPause;
    x.opts.iWorker = parameters.iWorker;
    x.opts.nWorkers = parameters.nWorkers;
    
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
    if ~x.opts.bProcessData
        reportProcess = 'run the initialization';
        if x.opts.bOnlyLoad
            reportProcess = 'load the dataset';
        end
    elseif x.opts.bProcessData
        reportProcess = 'run the processing pipeline';
    end
    if x.opts.bImportData
        reportImport = 'will run the import workflow and ';
    else
        reportImport = '';
    end
    % Print feedback
    fprintf('ExploreASL %swill %s...\n',reportImport,reportProcess);

end


%% -----------------------------------------------------------------------
function [x] = ExploreASL_Initialize_SubStructs(x)

    x.S = struct; % Statistics
    x.D = struct; % Directories
    x.P = struct; % Paths
    x.Q = struct; % Quality
    
    x.settings = struct;    % Workflow settings
    x.dataset = struct;     % Dataset related fields
    x.external = struct;    % Toolbox related fields (SPM, CAT, etc.)
    x.dir = struct;         % BIDS related directories (sourceStructure, studyPar, dataset_description, etc.)

end


