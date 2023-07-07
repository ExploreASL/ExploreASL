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
% This  workflow initializes ExploreASL:
% 
% 1. Admin (Validate Matlab version, parse input parameters, check which modules should be run)
% 2. Get main ExploreASL path
% 3. Add ExploreASL paths to Matlab and prepare x (sub)structs
% 4. Check DatasetRoot provided by user & provide feedback on which ExploreASL module is run
% 5. Define several general settings (mostly used for processing, PM: move these to processing initialization)
% 6. Print logo, settings, ExploreASL version etc
%
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:
%
% Calling externally:             x = ExploreASL_Initialize('/MyDisk/MyStudy');
% Debugging/initialization only:  x = ExploreASL_Initialize;
%
% __________________________________
% Copyright (c) 2015-2023 ExploreASL


    %% 1. Admin
    % Validate Matlab version, parse input parameters, check which modules should be run.

	if ~exist('x','var') || isempty(x)
		x = struct;
	end
    

    % Validate the Matlab version
	fprintf('%s\n', ['Matlab version detected: ' version]);
	fprintf('%s\n', 'Minimal Matlab version required: 9.6 (R2019a)');

	if verLessThan('matlab', '9.6')
		disp('<a href="https://nl.mathworks.com/support/requirements/previous-releases.html; ">Click here for the Matlab version overview</a>');
		error('Too old Matlab version for ExploreASL');
	end

    %% xASL_init_InputParsing Parse the input parameters of ExploreASL, 
    % meaning to divide and identify the input parameters, if they are provided in the correct format, 
    % define their defaults, and their meaning for ExploreASL.
    x = xASL_init_InputParsing(x,varargin{:});
	
    
    %% xASL_init_GetBooleansImportProcess Check which ExploreASL modules should be run, 
    % as per the user-provided booleans (e.g., bImport, bProcess)
    x = xASL_init_GetBooleansImportProcess(x);
    
    


    %% 2. Get main ExploreASL path

    x = xASL_init_GetMainExploreASLfolder(x);

    





    %% 3. Add ExploreASL paths to Matlab and prepare x (sub)structs


    %% xASL_init_AddDirsOfxASL  Add specific ExploreASL directories as Matlab paths, avoiding adding folders that are unused or that may even cause conflict.
    % If toolboxes that ExploreASL uses were already added to the paths, they are removed from the Matlab paths (such as spm, bids) to avoid conflicts.
    xASL_init_AddDirsOfxASL(x.opts.MyPath);

    
	% Initialize substructs
    % xASL_init_SubStructs Initialize the ExploreASL x structure substructs/fields.
    % Only fields which do not exist so far are added.
    x = xASL_init_SubStructs(x);
    



    %% 4. Check DatasetRoot provided by user & provide feedback on which ExploreASL module is run

    if x.opts.bImportData || x.opts.bLoadData || x.opts.bProcessData
        % xASL_init_checkDatasetRoot Check the provided ExploreASL input parameter "DatasetRoot".
        %  Is it a path? Does it exist? is it a directory? It should be the dataset root folder, 
        % but in previous versions of ExploreASL this used to be the path to the dataPar.json file.
        % This is managed here as well for backward compatibility (PM: backward compatibility can be removed)
        x = xASL_init_checkDatasetRoot(x);
    end
    


    % Provide feedback on which ExploreASL module is run
    % (PM: merge with xASL_init_printSettings below)
    xASL_init_BasicFeedback(x);
    



    %% 5. Define several general settings (mostly used for processing, PM: move these to processing initialization)

    % Go to ExploreASL folder
    cd(x.opts.MyPath);


    % These settings are data-independent
    % Define several data independent settings for processing 
    % (PM: move these to the processing initialization)
    x = xASL_init_DefineIndependentSettings(x);





    %% 6. Print logo, settings, ExploreASL version etc
    % xASL_init_PrintLogo Print the ExploreASL logo before running ExploreASL
    xASL_init_PrintLogo();


    % Print chosen settings
    % xASL_init_printSettings Print chosen settings, for:
    %  - Dataset root (which dataset root folder is selected)
    % - Import modules (which modules are run)
    % - Process modules (which modules are run)
    % - Pause before processing (yes/no)
    % - iWorker and nWorkers (for parallelization)
    xASL_init_printSettings(x);

    % Check permissions 
    % This part is skipped for now, this checked if the data and ExploreASL paths had sufficient permissions)
    % xASL_adm_CheckPermissions(x.opts.MyPath, false);


    % xASL_init_PrintVersion Print ExploreASL version to the screen, beta as well
    xASL_init_PrintVersion(x.Version);

	if ~isdeployed && usejava('desktop') % true if the Matlab GUI is loaded, false when in CLI with or without Java VM
        disp('<a href="https://exploreasl.github.io/Documentation; ">Click here for the ExploreASL manual</a>');
    else % text only
        fprintf('ExploreASL manual is available at https://exploreasl.github.io/Documentation\n');
    end
end




%% =======================================================================================================================
%% =======================================================================================================================
%% SUBFUNCTIONS start here 








%% 1. Admin ==============================================================================================================
%% =======================================================================================================================
%% =======================================================================================================================
function x = xASL_init_InputParsing(x, varargin)
    %% xASL_init_InputParsing Parse the input parameters of ExploreASL, 
    % meaning to divide and identify the input parameters, if they are provided in the correct format, 
    % define their defaults, and their meaning for ExploreASL.
    
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
    
    
    %% =======================================================================================================================
    %% =======================================================================================================================
    function x = xASL_init_GetBooleansImportProcess(x)
    % xASL_init_GetBooleansImportProcess Check which ExploreASL modules should be run, 
    % as per the user-provided booleans (e.g., bImport, bProcess)
    
    
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







%% 2. Get main ExploreASL path ===========================================================================================
%% =======================================================================================================================
%% =======================================================================================================================
function x = xASL_init_GetMainExploreASLfolder(x)
    % xASL_init_GetMainExploreASLfolder by
    % 1. Check if the current directory is the ExploreASL main directory
    % 2. Allow the user to input the ExploreASL root path manually
    % 3. In deployed mode set the ExploreASL directory in the ctf archive
    % 4. Go to ExploreASL directory
    
        %% 1. Check if the current directory is the ExploreASL main directory
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
    
        %% 2. Allow the user to input the ExploreASL root path manually
        MasterScriptPath = fullfile(x.opts.MyPath, 'ExploreASL.m');
    
        % Select the ExploreASL folder manually, if the script is not run in deployed mode
        if ~isdeployed
            if ~exist(MasterScriptPath,'file')
                pathstr = input('Provide foldername where ExploreASL is installed (format: \''PathExploreASL\''): ');
                if sum(pathstr==0) || ~exist(fullfile(pathstr,'ExploreASL.m'),'file'), return; end
                x.opts.MyPath = pathstr;
            end
        else
            %% 3. In deployed mode set the ExploreASL directory in the ctf archive
            
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
    
        %% 4. Go to ExploreASL directory
        cd(x.opts.MyPath);
    
    end




%% 3. Add ExploreASL paths ===============================================================================================
%% =======================================================================================================================
%% =======================================================================================================================
function xASL_init_AddDirsOfxASL(MyPath)
%% Add specific ExploreASL directories as Matlab paths, avoiding adding folders that are unused or that may even cause conflict.
% If toolboxes that ExploreASL uses were already added to the paths, they are removed from the Matlab paths (such as spm, bids) to avoid conflicts.

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
            fullfile('Modules', 'Module_Import'), ...
            fullfile('Modules', 'Module_Processing'), ...
            fullfile('Modules', 'Module_Processing', 'Module_Structural'), ...
            fullfile('Modules', 'Module_Processing', 'Module_ASL'), ...
            fullfile('Modules', 'Module_Processing', 'Module_Population'), ...
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









    
    

%% 4. Check DatasetRoot provided =========================================================================================
%% =======================================================================================================================
%% =======================================================================================================================


    
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
        % DESCRIPTION: Check the provided ExploreASL input parameter "DatasetRoot".
        %  Is it a path? Does it exist? is it a directory? It should be the dataset root folder, 
        % but in previous versions of ExploreASL this used to be the path to the dataPar.json file.
        % This is managed here as well for backward compatibility (PM: backward compatibility can be removed)
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
        
            % Track if a valid path was provided
            bValidPath = true;
        
            % Check if user correctly inserted a dataset root directory
            % Trying to fix incorrectly specified files
            [fPath, fFile, fExt] = fileparts(x.opts.DatasetRoot);
            if strcmp(fExt, '.json')
                warning([fFile fExt ' incorrectly provided as the dataset-root input. ExploreASL requires the path of the dataset root directory... Using ' fPath ' instead.']);
                x.opts.DatasetRoot = fPath;
            elseif ~isempty(fExt) && ~exist(x.opts.DatasetRoot, 'dir')
                % Files are not supported as dataset root directory. If the provided path does not exist as a directory (so it's potentially a file) and it has an extension, then we report that a file is incorrectly provided.
                % Note that a directory with '.' in the name can be incorrectly identified as having an 'extension' and confused with a file - therefore an extra check is presented here.
                warning([fFile fExt ' incorrectly provided as the dataset-root input. ExploreASL requires the path of the dataset root directory...']);
                bValidPath = false;
            end
        
            % Trying to fix incorrectly specified subdirectories
            % (note that this can depend on the above file check, e.g., in the case of /derivatives/ExploreASL/dataPar.json
            % so this should not be an elseif statement)
            if bValidPath
                [fPath, fSubFolder, fExt] = fileparts(x.opts.DatasetRoot);
                switch [fSubFolder fExt]
                    case 'sourcedata'
                        warning(['sourcedata directory provided as the dataset-root input. Using ' fPath ' instead.']);
                        x.opts.DatasetRoot = fPath;
                    case 'rawdata'
                        warning(['rawdata directory provided as the dataset-root input. Using ' fPath ' instead.']);
                        x.opts.DatasetRoot = fPath;
                    case 'derivatives'
                        warning(['derivatives directory provided as the dataset-root input. Using ' fPath ' instead.']);
                        x.opts.DatasetRoot = fPath;
                    case 'ExploreASL'
                        [fPath, fFolder] = fileparts(fPath);
                        if strcmp(fFolder, 'derivatives')
                            warning(['derivatives/ExploreASL directory provided as the dataset-root input. Using ' fPath ' instead.']);
                            x.opts.DatasetRoot = fileparts(fPath);
                        end
                end
            end
        
            % Check if the user provided any path, and if this path exists
            if isempty(x.opts.DatasetRoot)
                warning('Dataset root directory is not specified...');
                bValidPath = false;
            elseif ~exist(x.opts.DatasetRoot, 'dir')
                warning('Dataset root directory does not exist...');
                bValidPath = false;
            else
                % Define the other paths
                x = xASL_init_DetermineRequiredPaths(x);
            end    
        
            if ~bValidPath && (x.opts.bProcessData || x.opts.bImportData)
                % Give back a warning that the user tried to import or process but neither a correct dataset root nor a dataPar.json that exists was used
                warning('You are trying to import or process a dataset, but the input parameters are not correct. ExploreASL will only be initialized...');
                x.opts.bProcessData = 0;
                x.opts.bImportData = 0;
                x.opts.bLoadData = false;
                x.opts.bProcess = [0 0 0];
                x.opts.bImport = [0 0 0];
            end
            
            % Set defaults for "dir" fields
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



        %% =======================================================================================================================
        %% =======================================================================================================================
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
            
                % Regular expressions used below, explained:
                % ^ = expression must start with this
                % $ = expression must end with this
                % (?i) expression becomes case insensitive
                % \ = escaping
                % . = any character
                % * = previous repeated 0-Inf times
                % .* = any character repeated 0-Inf times
                % \. extension delimiter/period
                
                % BIDS DatasetRoot directory
                x.dir.DatasetRoot = x.opts.DatasetRoot;
                x.opts.dataParType = 'directory';
                
                % Set other basic BIDS directories
                x.dir.SourceData = fullfile(x.opts.DatasetRoot,'sourcedata');
                x.dir.RawData = fullfile(x.opts.DatasetRoot,'rawdata');
                x.dir.Derivatives = fullfile(x.opts.DatasetRoot,'derivatives');
                x.dir.xASLDerivatives = fullfile(x.dir.Derivatives,'ExploreASL');
                
                % Search for descriptive JSON files
                fileListSourceStructure = xASL_adm_GetFileList(x.dir.DatasetRoot, '(?i)^sourcestructure.*\.json$');
                fileListStudyPar = xASL_adm_GetFileList(x.dir.DatasetRoot, '(?i)^studypar.*\.json$');
                fileListDataDescription = xASL_adm_GetFileList(fullfile(x.dir.DatasetRoot, 'rawdata'), '(?i)^dataset_description\.json$');
                
                % First try the derivatives folder
                fileListDataPar = xASL_adm_GetFileList(fullfile(x.dir.DatasetRoot, 'derivatives', 'ExploreASL'), '(?i)^datapar.*\.json$');
                
                % Check for x.D.ROOT
                if exist(fullfile(x.dir.DatasetRoot, 'derivatives', 'ExploreASL'),'dir')
                    x.D.ROOT = fullfile(x.dir.DatasetRoot, 'derivatives', 'ExploreASL');
                end
                if isempty(fileListDataPar)
                    % Derivatives maybe does not exist already, we'll try study root
                    fileListDataPar = xASL_adm_GetFileList(x.dir.DatasetRoot, '(?i)^datapar.*\.json$');
                end
                
                % Assign fields
                if ~isempty(fileListSourceStructure)
                    if length(fileListSourceStructure) > 1
                        warning('Multiple sourcestructure*.jsons exist. Using the first: %s\n',fileListSourceStructure{1});
                    end
                    x.dir.sourceStructure = fileListSourceStructure{1};
                end
                if ~isempty(fileListStudyPar)
                    if length(fileListStudyPar) > 1
                        warning('Multiple studyPar*.jsons exist. Using the first: %s\n',fileListStudyPar{1});
                    end
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
            
            
            


        
%% =======================================================================================================================
%% =======================================================================================================================
function xASL_init_BasicFeedback(x)
    % xASL_init_BasicFeedback Provide feedback on which ExploreASL module is run
    % (PM: merge with the other feedback subfunction(s)?)
    
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



%% 5. Define several general settings ====================================================================================
%% =======================================================================================================================
%% =======================================================================================================================
    function [x] = xASL_init_DefineIndependentSettings(x)
        %xASL_init_DefineIndependentSettings Define ExploreASL environment parameters, independent of loaded data
        %
        % FORMAT: [x] = xASL_init_DefineIndependentSettings(x)
        %
        % INPUT:
        %   x       - ExploreASL x structure (STRUCT, REQUIRED)
        %
        % OUTPUT:
        %   x       - ExploreASL x structure (STRUCT)
        %
        % -----------------------------------------------------------------------------------------------------------------------------------------------------
        % DESCRIPTION: Define ExploreASL environment parameters, independent of loaded data.
        %
        % EXAMPLE:     This is part of the initialization workflow. Check out the usage there.
        %
        % -----------------------------------------------------------------------------------------------------------------------------------------------------
        % REFERENCES:  n/a
        %
        % Copyright 2015-2021 ExploreASL
        
        
        %% --------------------------------------------------------------------------------------------------------------------
        %% Default parameters
        x.settings.stopAfterErrors = Inf; % set to a high number (or Inf) when running a complete batch overnight
        x.settings.dryRun = false; % set to true to skip all real processing & just print all parameters
        x.settings.bOverwrite = true;
        
        if ~exist('groot','builtin')
            % before R2012b
            set(0,'DefaultTextInterpreter','none')
        else
            % R2012b+
            set(groot,'DefaultTextInterpreter','none')
        end
        
        % Get version
        try
            versionFile = dir(fullfile(x.opts.MyPath, 'VERSION*'));
            versionFile = versionFile.name;
            x.Version = versionFile(9:end);
        
            if isempty(x.Version)
                warning('Unknown ExploreASL version number');
                x.Version = '0.0.0 VERSION_unknown';
            end
        catch
            warning('Could not obtain ExploreASL version, version file missing');
            x.Version = '0.0.0 VERSION_File_Missing';
        end
        
        
        if ~isfield(x,'Q')
            x.Q = struct;
        end
        
        
        %% Atlases and templates
        x = xASL_init_MapsAndAtlases(x);
        
        %% Visualization
        x = xASL_init_VisualizationSettings(x);
        
        
        end    





%% 6. Print ==============================================================================================================
%% =======================================================================================================================
%% =======================================================================================================================
function xASL_init_PrintLogo()
    % xASL_init_PrintLogo Print the ExploreASL logo before running ExploreASL
    
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



%% =======================================================================================================================
%% =======================================================================================================================
function xASL_init_printSettings(x)
    %xASL_init_printSettings Print chosen settings
    %
    % FORMAT: xASL_init_printSettings(x)
    %
    % INPUT:
    %   x       - ExploreASL x structure (STRUCT, REQUIRED)
    %
    % OUTPUT:
    %   n/a
    %
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % DESCRIPTION: Print chosen settings, for:
    %  - Dataset root (which dataset root folder is selected)
    % - Import modules (which modules are run)
    % - Process modules (which modules are run)
    % - Pause before processing (yes/no)
    % - iWorker and nWorkers (for parallelization)
    %
    % EXAMPLE:     This is part of the initialization workflow. Check out the usage there.
    %
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % REFERENCES:  n/a
    %
    % Copyright 2015-2022 ExploreASL
    
        %% Printing
        xASL_adm_BreakString('ExploreASL Settings');
        % Dataset root
        if length(x.opts.DatasetRoot)>70
            fprintf('Dataset Root        ...%s\n', x.opts.DatasetRoot(end-70:end));
        else
            fprintf('Dataset Root        %s\n', x.opts.DatasetRoot);
        end
    
        % Import modules
        textPrint = 'Import & Defacing   ';
        if x.opts.bImport(1)==1
            textPrint = [textPrint 'DCM2NII '];
        end
        if x.opts.bImport(2)==1
            textPrint = [textPrint 'NII2BIDS '];
        end
        if x.opts.bImport(3)==1
            textPrint = [textPrint 'DEFACE '];
        end
        fprintf([textPrint '\n']);
    
        % Process modules
        textPrint = 'Processing          ';
        if x.opts.bProcess(1)==1
            textPrint = [textPrint 'STRUCTURAL '];
        end
        if x.opts.bProcess(2)==1
            textPrint = [textPrint 'ASL '];
        end
        if x.opts.bProcess(3)==1
            textPrint = [textPrint 'POPULATION '];
        end
        fprintf([textPrint '\n']);
    
        % Pause before processing
        if x.opts.bPause==1
            fprintf('bPause              %s\n', 'True');
        else
            fprintf('bPause              %s\n', 'False');
        end
    
        % Worker numbers
        fprintf('iWorker             %d\n', x.opts.iWorker);
        fprintf('nWorkers            %d\n', x.opts.nWorkers);
        xASL_adm_BreakString('',[],[],1);
    
    end


%% =======================================================================================================================
%% =======================================================================================================================
function xASL_init_PrintVersion(vExploreASL)
    % Print ExploreASL version to the screen, beta as well
    
        % For beta versions we print the text in red with an additional comment, to make sure users are aware of it
        if isempty(regexpi(vExploreASL,'BETA'))
            fprintf('ExploreASL v%s initialized ... \n', vExploreASL);
        else
            vExploreASL = vExploreASL(1:regexpi(vExploreASL,'_BETA')-1);
            fprintf(2,'ExploreASL v%s (beta) initialized... \n', vExploreASL);
        end
    
    end