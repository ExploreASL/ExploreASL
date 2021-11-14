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
% 5. General settings
% - Check if the dataPar.json exists and if it can be loaded or not
% - Initialize the data independent & dependent settings
%
% 6. Print logo & settings, check permissions
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
    x.settings.SelectParFile = false; % Fallback
    if x.opts.bProcessData
        % Checkout the "Proceed with Initialization" section
        if (isempty(x.opts.DatasetRoot) || (~exist(x.opts.DatasetRoot,'file') && ~exist(x.opts.DatasetRoot,'dir')))
            x.settings.SelectParFile = true; % If the DatasetRoot is either empty OR the file does not exist, we have to select it later on (if processing is turned on)
        end
    end
    
    
    %% 2. Get ExploreASL path

    % Check if the current directory is the ExploreASL directory
    CurrCD = pwd;
    if exist(fullfile(CurrCD, 'ExploreASL_Master.m'), 'file')
        x.opts.MyPath = CurrCD;
    end

    % Check whether MyPath is correct, otherwise obtain correct folder
    if ~isfield(x.opts, 'MyPath')
        % Check if we can get the path from the ExploreASL_Initialize path
        initializePath = fileparts(mfilename('fullpath'));
        if ~isempty(regexp(initializePath,'ExploreASL$', 'once'))
            x.opts.MyPath = initializePath;
        else
            x.opts.MyPath = '/DummyPath';
        end
    end

    % Get the master script path
    MasterScriptPath = fullfile(x.opts.MyPath, 'ExploreASL_Master.m');

    % Select the ExploreASL folder manually, if the script is not run in deployed mode
    if ~isdeployed
        if ~exist(MasterScriptPath,'file')
            pathstr = input('Provide foldername where ExploreASL is installed (format: \''PathExploreASL\''): ');
            if sum(pathstr==0) || ~exist(fullfile(pathstr,'ExploreASL_Master.m'),'file'), return; end
            x.opts.MyPath = pathstr;
        end
    else
        % In deployed mode set the ExploreASL directory in the ctf archive
        [files,~] = spm_select('FPListRec',ctfroot,'ExploreASL_Master*'); % Find the path of the master files within the ctf archive
        curPathCTF = fileparts(files(1,:)); % Get the path
        x.opts.MyPath = fullfile(curPathCTF); % curPathCTF = ExploreASL path

        BreakString = '==============================================================================================\n';
        fprintf(BreakString);
        fprintf('ctfroot:  %s\n', ctfroot);
        fprintf('x.opts.MyPath: %s\n', x.opts.MyPath);
        fprintf(BreakString);

    end

    % Go to ExploreASL folder
    cd(x.opts.MyPath);


    %% 3. Add ExploreASL paths
    if ~isdeployed
        addExploreASLDirectory(x.opts.MyPath);
    end
    
    
    %% 4. Check DatasetRoot
    [x] = xASL_init_checkDatasetRoot(x);
    
    % Give some feedback
    ExploreASL_Initialize_basicFeedback(x);
    
    %% 5. General settings

    % Go to ExploreASL folder
    cd(x.opts.MyPath);

    % These settings are data-independent
    x = xASL_init_DefineIndependentSettings(x);

    %% 6. Print logo & settings, check permissions
    xASL_init_PrintLogo(x);

    % Print chosen settings
    xASL_init_printSettings(x);

    % Check permissions
    % xASL_adm_CheckPermissions(x.opts.MyPath, false);

    %% 7. Data-specific initialization
    xASL_init_PrintVersion(x.Version);

end


%% ==================================================================================
function xASL_init_PrintVersion(vExploreASL)

    % For beta versions we print the text in red with an additional comment, to make sure users are aware of it
    if isempty(regexpi(vExploreASL,'BETA'))
        fprintf('ExploreASL v%s initialized ... \n', vExploreASL);
    else
        vExploreASL = vExploreASL(1:regexpi(vExploreASL,'_BETA')-1);
        fprintf(2,'ExploreASL v%s (beta) initialized... \n', vExploreASL);
    end

end


%% ==================================================================================
function xASL_init_PrintLogo(x)

    % Add design fields
    logoString = sprintf([...
    '[\b' ...
    ' ________                      __                                 ______    ______   __          \n'...
    '/        |                    /  |                               /      \\  /      \\ /  |       \n'...
    '########/  __    __   ______  ## |  ______    ______    ______  /######  |/######  |## |         \n'...
    '## |__    /  \\  /  | /      \\ ## | /      \\  /      \\  /      \\ ## |__## |## \\__##/ ## |   \n'...
    '##    |   ##  \\/##/ /######  |## |/######  |/######  |/######  |##    ## |##      \\ ## |       \n'...
    '#####/     ##  ##<  ## |  ## |## |## |  ## |## |  ##/ ##    ## |######## | ######  |## |         \n'...
    '## |_____  /####  \\ ## |__## |## |## \\__## |## |      ########/ ## |  ## |/  \\__## |## |_____ \n'...
    '##       |/##/ ##  |##    ##/ ## |##    ##/ ## |      ##       |## |  ## |##    ##/ ##       |   \n'...
    '########/ ##/   ##/ #######/  ##/  ######/  ##/        #######/ ##/   ##/  ######/  ########/    \n'...
    '                    ## |                                                                         \n'...
    '                    ## |                                                                         \n'...
    '                    ##/  ]\b\n\n']);

    xASL_adm_BreakString('');
    fprintf(logoString);

end


%% -----------------------------------------------------------------------
%% Add ExploreASL Directory
function addExploreASLDirectory(MyPath)

    % Define paths (should be equal when loading data or only initializing)
    addpath(MyPath); % ExploreASL
    subfoldersToAdd = {...
        'External', 'Functions', 'mex', 'Modules', 'Testing',...
        fullfile('Development', 'Outdated'),...
        fullfile('Modules', 'SubModule_Import'), ...
        fullfile('Modules', 'SubModule_Structural'), ...
        fullfile('Modules', 'SubModule_ASL'), ...
        fullfile('Modules', 'SubModule_Population'), ...
        fullfile('External','isnear'), ...
        fullfile('External','DCMTK'), ...
        fullfile('Testing', 'UnitTests'), ...
        fullfile('Testing', 'Functions'), ...
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
    if length(parameters.ImportModules)==1
        % If a single value is given ...
        % ... then turn on/off all the import modules ...
        parameters.ImportModules(1:4) = logical(parameters.ImportModules(1));
        % ... besides defacing, which can only be run using a 1x4 vector ...
        parameters.ImportModules(3) = false;
    elseif length(parameters.ImportModules)<4
        % Convert to a row vector
        parameters.ImportModules = parameters.ImportModules(:)';
        % Fill in the missing fields with zeros
        parameters.ImportModules(length(parameters.ImportModules)+1:4) = 0;
        % Issue a warning
        warning('Incorrect length of the ImportModules parameter, missing submodules set to zero: %s\n',num2str(parameters.ImportModules));
    end
    if length(parameters.ProcessModules)==1
        % If a single value is given, then copy it to all submodules
        parameters.ProcessModules(1:3) = parameters.ProcessModules(1);
    elseif length(parameters.ProcessModules)<3
        % Convert to a row vector
        parameters.ProcessModules = parameters.ProcessModules(:)';
        % Fill in the missing fields with zeros
        parameters.ProcessModules(length(parameters.ProcessModules)+1:3) = 0;
        % Issue a warning
        warning('Incorrect length of the ProcessModules parameter, missing submodules set to zero: %s\n',num2str(parameters.ProcessModules));
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

    % Check if data is being imported
    if sum(x.opts.ImportModules)>0
        x.opts.bImportData = 1;
    else
        x.opts.bImportData = 0;
    end
    
    % Check if data is being processed
    if sum(x.opts.ProcessModules)>0
        x.opts.bProcessData = 1;
    else
        x.opts.bProcessData = 0;
    end
    
    % We only load data when a DatasetRoot is provided
    if ~isempty(x.opts.DatasetRoot)
        x.opts.bLoadData = true;
    else
        x.opts.bLoadData = false;
    end

end


%% -----------------------------------------------------------------------
% Give some feedback
function ExploreASL_Initialize_basicFeedback(x)

    % Report string
    reportProcess = '';
    if ~x.opts.bProcessData
        reportProcess = 'run the initialization';
        if x.opts.bLoadData
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
    
    % Statistics, directories, paths, and sequence related fields
    if ~isfield(x,'S'),                     x.S = struct;                   end
    if ~isfield(x,'D'),                     x.D = struct;                   end
    if ~isfield(x,'P'),                     x.P = struct;                   end
    if ~isfield(x,'Q'),                     x.Q = struct;                   end
    
    % Module subfields (import, structural, asl, & population)
    if ~isfield(x,'modules'),               x.modules = struct;             end
    if ~isfield(x.modules,'import'),        x.modules.import = struct;      end
    if ~isfield(x.modules,'structural'),    x.modules.structural = struct;  end
    if ~isfield(x.modules,'asl'),           x.modules.asl = struct;         end
    if ~isfield(x.modules,'population'),    x.modules.population = struct;  end
    
    % Dataset related fields, workflow settings, toolbox/external (SPM, CAT, FSL, etc.) fields, general directories
    if ~isfield(x,'dataset'),               x.dataset = struct;             end
    if ~isfield(x,'settings'),              x.settings = struct;            end
    if ~isfield(x,'external'),              x.external = struct;            end
    if ~isfield(x,'dir'),                   x.dir = struct;                 end

end


