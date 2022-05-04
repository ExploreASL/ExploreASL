function [x] = ExploreASL_Initialize(varargin)
%ExploreASL_Initialize Initializes ExploreASL
%
% FORMAT: 
%   x = ExploreASL_Initialize([DatasetRoot, bImport, bProcess, bPause, iWorker, nWorkers])
%
% INPUT:
%   varargin    - This script accepts the same arguments as ExploreASL. Check out the definitions there.
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
% - Input parsing (xASL_init_InputParsing, xASL_init_convertParsedInput, xASL_init_storeParsedInput)
% - Initialize substructs of the ExploreASL x structure (xASL_init_SubStructures)
% - Check input parameters and determine pipeline related booleans (xASL_init_GetBooleansImportProcess)
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
% Calling externally:             x = ExploreASL_Initialize('/MyDisk/MyStudy');
% Debugging/initialization only:  x = ExploreASL_Initialize;
%
% __________________________________
% Copyright (c) 2015-2022 ExploreASL


    %% 1. Admin

	if ~exist('x','var') || isempty(x)
		x = struct;
	end
        
    % Parse input parameters
    x = xASL_init_InputParsing(x,varargin{:});
	    
    % Check if the ExploreASL pipeline should be run or not
    x = xASL_init_GetBooleansImportProcess(x);
    
    
    %% 2. Get ExploreASL path

    % Determine the MyPath of ExploreASL
    x = xASL_init_GetMyPath(x);

    % Get the master script path
    x = xASL_init_GetMasterScript(x);

    % Go to ExploreASL directory
    cd(x.opts.MyPath);

    %% 3. Add ExploreASL paths
    xASL_init_AddDirsOfxASL(x.opts.MyPath);
    
	% Initialize substructs
	x = xASL_init_SubStructs(x);
    
    %% 4. Check DatasetRoot
    x = xASL_init_checkDatasetRoot(x);
    
    % Give some feedback
    xASL_init_BasicFeedback(x);
    
    %% 5. General settings

    % Go to ExploreASL folder
    cd(x.opts.MyPath);

    % These settings are data-independent
    x = xASL_init_DefineIndependentSettings(x);

    %% 6. Print logo & settings, check permissions
    xASL_init_PrintLogo();

    % Print chosen settings
    xASL_init_printSettings(x);

    % Check permissions
    % xASL_adm_CheckPermissions(x.opts.MyPath, false);

    %% 7. Data-specific initialization
    xASL_init_PrintVersion(x.Version);

	if ~isdeployed && usejava('desktop') % true if the Matlab GUI is loaded, false when in CLI with or without Java VM
        disp('<a href="https://exploreasl.github.io/Documentation; ">Click here for the ExploreASL manual</a>');
    else % text only
        fprintf('ExploreASL manual is available at https://exploreasl.github.io/Documentation\n');
    end
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
function xASL_init_PrintLogo()

    % Add design fields
    logoString = [...
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
    '                    ##/                                                                          \n'];

    xASL_adm_BreakString('',[],[],1);
    fprintf(logoString);

end


%% -----------------------------------------------------------------------
%% Add ExploreASL Directory
function xASL_init_AddDirsOfxASL(MyPath)

    if ~isdeployed
        
        % First remove existing toolbox initializations
        % This could contain other toolbox versions and create conflicts
        currentPathList = path;
        if ischar(currentPathList)
            % Convert to cell
            indicesAre = [0 strfind(currentPathList, pathsep)];
            % Add to pathList
            for iIndex=1:numel(indicesAre)-1
                pathList{iIndex, 1} = currentPathList(indicesAre(iIndex)+1:indicesAre(iIndex+1)-1);
            end
        else
            pathList = currentPathList;
        end

        % Iterate over paths
        for iPath=1:numel(pathList)
            if ~isempty(regexpi(pathList{iPath}, [filesep '(spm|bids-matlab)'])) % toolboxes can be added here
                % Here we want search for toolboxes that ExploreASL uses, but
                % that are in another path (e.g., within the Matlab toolboxes)

                if isempty(regexp(pathList{iPath}, fullfile('ExploreASL', 'External'), 'once'))
                    % If this path is not an ExploreASL-contained toolbox
                    rmpath(pathList{iPath});
                    fprintf(2,'Warning: Removed Matlab path to avoid conflicts: %s\n',pathList{iPath});
                end
            end
        end
        if numel(pathList)>1
            fprintf('\n');
        end

        % Define paths (should be equal when loading data or only initializing)
        addpath(MyPath); % ExploreASL
        subfoldersToAdd = {...
            'External', 'Modules', 'Testing', 'Functions',...
            fullfile('Functions', 'Administration'), ...
            fullfile('Functions', 'BIDS'), ...
		    fullfile('Functions', 'Development'), ...			
            fullfile('Functions', 'FSL'), ...
            fullfile('Functions', 'ImageProcessing'), ...
            fullfile('Functions', 'Import'), ...
            fullfile('Functions', 'Initialization'), ...
            fullfile('Functions', 'InputOutput'), ...
            fullfile('Functions', 'mex'), ...
            fullfile('Functions', 'QualityControl'), ...
            fullfile('Functions', 'Quantification'), ...
            fullfile('Functions', 'SPM'), ...
            fullfile('Functions', 'Statistics'), ...
            fullfile('Functions', 'Visualization'), ...
            fullfile('Modules', 'SubModule_Import'), ...
            fullfile('Modules', 'SubModule_Structural'), ...
            fullfile('Modules', 'SubModule_ASL'), ...
            fullfile('Modules', 'SubModule_Population'), ...
            fullfile('External','isnear'), ...
            fullfile('External','DCMTK'), ...
            fullfile('Testing', 'UnitTests'), ...
            fullfile('Testing', 'Functions'), ...
            fullfile('External','ExploreQC'), ...
            fullfile('External','bids-matlab'), ...
            fullfile('External','SPMmodified'), ...
            fullfile('External','SPMmodified','matlabbatch'),...
            fullfile('External','SPMmodified','xASL'),...
            fullfile('External','SPMmodified','toolbox','cat12'), ...
            fullfile('External','SPMmodified','toolbox','LST'), ...
            fullfile('External','SPMmodified','toolbox','OldNorm')};
            % ...
            % genpath(fullfile('External','bids-matlab'))};
        % Iterate over subfolders which should be added
        for ii=1:length(subfoldersToAdd)
            addpath(fullfile(MyPath,subfoldersToAdd{ii}));
        end
    end
    
end

%% -----------------------------------------------------------------------
%% Define input parser
function x = xASL_init_InputParsing(x, varargin)

    % This subfunction receives the input arguments and parses them.
    % The allowed input arguments are: 
	%     DatasetRoot - a path to a directory or an empty variable
	%     bImport     - a boolean in the form of a char or numerical or logical array
	%     bProcess    - a boolean in the form of a char or numerical or logical array
	%     bPause      - a boolean in the form of a char or numerical or logical value
	%     iWorker     - a boolean in the form of a char or numerical value
	%     nWorkers    - a boolean in the form of a char or numerical value
	% The outputs are written into the x.opts structure

	% DatasetRoot
	if length(varargin)<1 || isempty(varargin{1})
		x.opts.DatasetRoot = '';
	elseif ~ischar(varargin{1})
		error('Parameter DatasetRoot does not have a valid format');
	else
		x.opts.DatasetRoot = varargin{1};
	end

	% bImport
	if length(varargin)<2 || isempty(varargin{2})
		x.opts.bImport = [0 0 0];
	elseif ~( ischar(varargin{2}) || isnumeric(varargin{2}) || islogical(varargin{2}) )
		error('Parameter bImport does not have a valid format');
	else
		x.opts.bImport = varargin{2};
		if ischar(x.opts.bImport)
			x.opts.bImport = str2num(x.opts.bImport);
		end
	end
	
	% bProcess
	if length(varargin)<3 || isempty(varargin{3})
		x.opts.bProcess = [0 0 0];
	elseif ~( ischar(varargin{3}) || isnumeric(varargin{3}) || islogical(varargin{3}) )
		error('Parameter bProcess does not have a valid format');
	else
		x.opts.bProcess = varargin{3};
		if ischar(x.opts.bProcess)
			x.opts.bProcess = str2num(x.opts.bProcess);
		end
	end
	
	% bPause
	if length(varargin)<4 || isempty(varargin{4})
		x.opts.bPause = 0;
	elseif ~( ischar(varargin{4}) || isnumeric(varargin{4}) || islogical(varargin{4}) )
		error('Parameter bPause does not have a valid format');
	else
		x.opts.bPause = varargin{4};
		if ischar(x.opts.bPause)
			x.opts.bPause = str2num(x.opts.bPause);
		end
	end
	
	% iWorker
	if length(varargin)<5 || isempty(varargin{5})
		x.opts.iWorker = 1;
	elseif ~( ischar(varargin{5}) || isnumeric(varargin{5}) )
		error('Parameter iWorker does not have a valid format');
	else
		x.opts.iWorker = varargin{5};
		if ischar(x.opts.iWorker)
			x.opts.iWorker = str2num(x.opts.iWorker);
		end
	end
	
	% nWorkers
	if length(varargin)<6 || isempty(varargin{6})
		x.opts.nWorkers = 1;
	elseif ~( ischar(varargin{6}) || isnumeric(varargin{6}) )
		error('Parameter nWorkers does not have a valid format');
	else
		x.opts.nWorkers = varargin{6};
		if ischar(x.opts.nWorkers)
			x.opts.nWorkers = str2num(x.opts.nWorkers);
		end
	end
    
    % Check length of arrays (single digit input)
    if length(x.opts.bImport)==1
        % If a single value is given then turn on/off all the import modules, but not defacing
        x.opts.bImport(1:2) = logical(x.opts.bImport(1));
		x.opts.bImport(3) = false;
    elseif length(x.opts.bImport)==4
        % Old version
        warning('You are using the outdated format with 4 import modules, additional elements are skipped...');
        x.opts.bImport = logical(x.opts.bImport(1:3));
    elseif length(x.opts.bImport)<3
        % Convert to a row vector
        x.opts.bImport = logical(x.opts.bImport(:)');
        % Issue a warning
        warning('Incorrect number of import modules (%s), missing sub-modules set to false...', num2str(length(x.opts.bImport)));
        % Fill in the missing fields with zeros
        x.opts.bImport(length(x.opts.bImport)+1:3) = false;
    elseif length(x.opts.bImport)>3
        % Skip additional elements
        warning('Incorrect number of import modules (%s), additional elements are skipped...', num2str(length(x.opts.bImport)));
        x.opts.bImport = logical(x.opts.bImport(1:3));
	end
	
    if length(x.opts.bProcess)==1
        % If a single value is given, then copy it to all submodules
        x.opts.bProcess(1:3) = logical(x.opts.bProcess(1));
    elseif length(x.opts.bProcess)<3
        % Convert to a row vector
        x.opts.bProcess = logical(x.opts.bProcess(:)');
        % Issue a warning
        warning('Incorrect number of processing modules (%s), missing sub-modules set to zero...', num2str(length(x.opts.bProcess)));
        % Fill in the missing fields with false
        x.opts.bProcess(length(x.opts.bProcess)+1:3) = false;
    elseif length(x.opts.bProcess)>3
        % Skip additional elements
        warning('Incorrect number of processing modules (%s), additional elements are skipped...', num2str(length(x.opts.bProcess)));
        x.opts.bProcess = logical(x.opts.bProcess(1:3));
    end
    
    % Make it impossible to set bPause to true in deployed mode
    if isdeployed
        x.opts.bPause = false;
    end
    
    % Check nWorkers for the parallelization of ExploreASL
    if numel(x.opts.nWorkers)>1
        warning('nWorkers is supposed to be a scalar, not a vector, resetting to the first element...');
        x.opts.nWorkers = x.opts.nWorkers(1);
    end
    
end


%% -----------------------------------------------------------------------
% Check if the ExploreASL pipeline should be run or not
function x = xASL_init_GetBooleansImportProcess(x)

    % On default we assume that we cannot load the data
    x.opts.bLoadableData = false;
    
    % Field to check if the data was loaded or not
    x.opts.bDataLoaded = false;

    % Check if data is being imported
    if sum(x.opts.bImport)>0
        x.opts.bImportData = 1;
    else
        x.opts.bImportData = 0;
    end

    % Check if data is being processed
    if sum(x.opts.bProcess)>0
        x.opts.bProcessData = 1;
    else
        x.opts.bProcessData = 0;
    end
    
    % We only load data when a DatasetRoot is provided
    if ~isempty(x.opts.DatasetRoot)
        % We only load the data if bImportData and bProcessData are false or only bProcessData is true
        if (~x.opts.bImportData && ~x.opts.bProcessData) || x.opts.bProcessData
            x.opts.bLoadData = true;
        end
    else
        x.opts.bLoadData = false;
    end
    
    % We need to check if import/defacing and processing were run separately
    if x.opts.bImportData && ~x.opts.bProcessData
        x.opts.bSkipBIDS2Legacy = true;
        x.opts.bLoadData = false;
    else
        x.opts.bSkipBIDS2Legacy = false;
    end

    % Warn the user if processing data but not loading a dataset
    if x.opts.bProcessData && ~x.opts.bLoadData
        warning('Data loading was disabled or incorrect dataset provided!');
        fprintf('Forcing data loading for processing now, this may go wrong\n');
        x.opts.bLoadData = true;
    end

end


%% -----------------------------------------------------------------------
% Give some feedback
function xASL_init_BasicFeedback(x)

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
function x = xASL_init_DeployedModeHandling(x)

    % Find the path of the master files within the ctf archive
    [files,~] = spm_select('FPListRec',ctfroot,'ExploreASL*');
    
    % Get the path
    curPathCTF = fileparts(files(1,:));
    
    % curPathCTF = ExploreASL path
    x.opts.MyPath = fullfile(curPathCTF);

    % User feedback
    BreakString = '==============================================================================================\n';
    fprintf(BreakString);
    fprintf('ctfroot:       %s\n', ctfroot);
    fprintf('x.opts.MyPath: %s\n', x.opts.MyPath);
    fprintf(BreakString);
    
end


%% -----------------------------------------------------------------------
function x = xASL_init_GetMyPath(x)

    % Check if the current directory is the ExploreASL directory
    CurrCD = pwd;
    if exist(fullfile(CurrCD, 'ExploreASL.m'), 'file')
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

end


%% -----------------------------------------------------------------------
function x = xASL_init_GetMasterScript(x)

    MasterScriptPath = fullfile(x.opts.MyPath, 'ExploreASL.m');

    % Select the ExploreASL folder manually, if the script is not run in deployed mode
    if ~isdeployed
        if ~exist(MasterScriptPath,'file')
            pathstr = input('Provide foldername where ExploreASL is installed (format: \''PathExploreASL\''): ');
            if sum(pathstr==0) || ~exist(fullfile(pathstr,'ExploreASL.m'),'file'), return; end
            x.opts.MyPath = pathstr;
        end
    else
        % In deployed mode set the ExploreASL directory in the ctf archive
        x = xASL_init_DeployedModeHandling(x);
    end

end
