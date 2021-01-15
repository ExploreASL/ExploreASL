function [x] = ExploreASL_Initialize(DataParPath, ProcessData, iWorker, nWorkers)
%ExploreASL_Initialize Initializes ExploreASL
%
% FORMAT: [x] = ExploreASL_Initialize([DataParPath, ProcessData, iWorker, nWorkers])
%
% INPUT:
%   DataParPath - path to data parameter file (OPTIONAL, required when ProcessData=true, will then be prompted if omitted)
%   ProcessData - 0 = only initialize ExploreASL functionality (e.g. path management etc, REQUIRED, will be prompted if omitted)
%               - 1 = initialize ExploreASL functionality, load data & start processing pipeline, 
%               - 2 = initialize ExploreASL functionality, load data but no processing
%               - (OPTIONAL, default=prompt the user)
%   iWorker     - allows parallelization when called externally. iWorker defines which of the parallel ExploreASL calls we are (OPTIONAL, DEFAULT=1)
%   nWorkers    - allows parallelization when called externally. nWorkers defines how many ExploreASL calls are made in parallel (OPTIONAL, DEFAULT=1)
%
% OUTPUT:
%   x           - struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This initialization wrapper initializes ExploreASL: managing paths, deployment, etc.
% Using the following initialization functions functions:
% xASL_init_DefinePaths            - manages folders for ExploreASL, and sets and creates study/data-folders
% xASL_init_Toolboxes              - initialization third-party toolboxes, e.g. SPM, dip_image (soon to be removed)
% xASL_init_VisualizationSettings  - defines visualization settings for
%                                    visual QC figure printing (type help xASL_init_VisualizationSettings for more information)
% xASL_init_DefineSets             - Define study subjects/parameters for this pipeline run
% xASL_init_PrintCheckSettings     - prints summarized data parameters and warnings
% xASL_init_FileSystem             - dirty initialization of common filenames used throughout the pipeline
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE for GUI: ExploreASL_Initialize
% EXAMPLE for calling externally: ExploreASL_Initialize('//MyDisk/MyStudy/DataPar.m', true, 1, 1);
% EXAMPLE for debugging/initialization only: [x] = ExploreASL_Initialize('',0);
% __________________________________
% Copyright 2015-2020 ExploreASL



%% Admin

% Defaults
if nargin<4 || isempty(nWorkers)
    nWorkers = 1;
end
if nargin<3 || isempty(iWorker)
    iWorker = 1;
end

x.nWorkers = 1;
x.iWorker = 1;

if nargin>3 && ~isempty(nWorkers)
    if isnumeric(nWorkers)
        x.nWorkers = nWorkers;
    else
        warning('argin "nWorkers" should be numerical, using default==1');
    end
end
if nargin>2 && ~isempty(iWorker)
    if isnumeric(iWorker)
        x.iWorker = iWorker;
    else
        warning('argin "iWorker" should be numerical, using default==1');
    end
end

if nargin<1 || isempty(DataParPath)
    SelectParFile = true;
elseif exist(DataParPath,'dir')
    SelectParFile = true;
    warning('Data parameter file was inputted as folder, but should point to a file instead');
else
    SelectParFile = false;
end

x.S = struct;

if nargin>1 && ~isempty(ProcessData)
    if ProcessData==0 || ProcessData==1  || ProcessData==2
        x.ProcessData = ProcessData; % skip GUI
    else
        warning('Invalid input argument "ProcessData", should be 0, 1, or 2, resetting');
    end
end


%% Question dialogue - initialize only or run data
 if ~usejava('desktop') || ~usejava('jvm') || ~feature('ShowFigureWindows')
     bUseGUI = false;
 else
     bUseGUI = true;
 end

if ~isfield(x,'ProcessData')

    InitChoice = 'Process dataset';
    if ~exist('DataParPath','var')
        ParFileExists = false;
    elseif isempty(DataParPath)
        ParFileExists = false;
    else
        ParFileExists = true;
    end

    if ~ParFileExists
        if bUseGUI
            InitChoice = questdlg('Would you like to initialize ExploreASL functionality, load a dataset and/or start the processing pipeline?', ...
                'Start up ExploreASL', 'Process dataset', 'Only initialize ExploreASL functionality', 'Load dataset only', 'Only initialize ExploreASL functionality');
            if isempty(InitChoice)
                InitChoice = 'Cancel';
            end
        else
            fprintf('Would you like to load a dataset or only initialize ExploreASL (set paths etc)?\n');
            cliChoice = input('Please press [1] for "Process dataset", [2] for "Only initialize ExploreASL functionality", [3] for "Load dataset only", or [4] to cancel: ');

            switch cliChoice
                case 1
                    InitChoice = 'Process dataset';
                case 2
                    InitChoice = 'Only initialize ExploreASL functionality';
                case 3
                    InitChoice = 'Load dataset only';
                case 4
                    InitChoice = 'Cancel';
                otherwise
                    fprintf('Wrong choice, please choose 1 2, 3, or 4: ');
                    return;
            end
        end
    end
    % Handle response
    switch InitChoice
        case 'Process dataset'
            fprintf('%s\n','Loading & processing dataset')
            x.ProcessData = 1;
        case 'Only initialize ExploreASL functionality'
            x.ProcessData = 0;
        case 'Load dataset only'
            fprintf('%s\n','Initializing ExploreASL functionality & loading dataset');
            x.ProcessData = 2;
        case 'Cancel'
            fprintf('%s\n','Exiting ExploreASL');
            x.ProcessData = 0;
            return;
        otherwise
            x.ProcessData = 0;
            fprintf('Unknown option, exiting\n');
            return;
    end
end

if x.ProcessData==2 && nargout==0
    warning('Data loading requested but no output structure defined');
    fprintf('%s\n', 'Try adding "x = " to the command to load data into the x structure');
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

MasterScriptPath = fullfile(x.MyPath, 'ExploreASL_Master.m');

if ~isdeployed
    if ~exist(MasterScriptPath,'file')
        if bUseGUI
            if ismac % somehow macOS doesnt show the title of the question dialogue
                fprintf('Select folder where ExploreASL is installed\n');
            end
            pathstr = uigetdir(CurrCD, 'Select folder where ExploreASL is installed');
        else
            pathstr = input('Provide foldername where ExploreASL is installed, including ": ');
        end
        if sum(pathstr==0) || ~exist(fullfile(pathstr,'ExploreASL_Master.m'),'file')
            return
        end
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
% Define paths (should be equal when loading data or only initializing
% ExploreASL
    addpath(x.MyPath);

	%fullfile('Development','dicomtools'), ...
    subfolders_to_add = {'Functions', 'mex', 'Modules', ...
                            fullfile('Modules', 'SubModule_Structural'), ...
                            fullfile('Modules', 'SubModule_ASL'), ...
                            fullfile('Modules', 'SubModule_Population'), ...
                            'Development', 'External', ...
                            fullfile('External','isnear'), ...
                            fullfile('External','DCMTK'), ...
                            fullfile('External','ExploreQC'), ...
							fullfile('External','SPMmodified'), ...
                            fullfile('External','SPMmodified','matlabbatch'),...                            
                            fullfile('External','SPMmodified','xASL'),...
                            fullfile('External','SPMmodified','toolbox','cat12'), ...
							fullfile('External','SPMmodified','toolbox','LST'), ...
							fullfile('External','SPMmodified','toolbox','OldNorm')};

    for ii=1:length(subfolders_to_add)
        addpath(fullfile(x.MyPath,subfolders_to_add{ii}));
    end
end


if x.ProcessData>0
    x = xASL_init_LoadDataParameterFile(x, DataParPath, SelectParFile, bUseGUI);
end


%% Initialize general parameters
x = xASL_init_DefineIndependentSettings(x); % these settings are data-independent

x = xASL_init_DefineDataDependentSettings(x); % these settings depend on the data (e.g. which template to use)


%% --------------------------------------------------------------------------------------------------------------------
%% Print logo
fprintf('%s\n',['--> Initializing ExploreASL v' x.Version '...']);

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

fprintf('\n\n');
fprintf([BreakString LogoString BreakString]);
fprintf('\n\n');



%% Check permissions
%xASL_adm_CheckPermissions(x.MyPath, false);




%% Data-specific initialization
fprintf('%s\n','ExploreASL initialized <--');

if x.ProcessData>0
    if isempty(x.D.ROOT)
        error('No root/analysis/study folder defined');
    end

    % Fix a relative path
    if strcmp(x.D.ROOT(1), '.')
        cd(x.D.ROOT);
        x.D.ROOT = pwd;
    end

    %% Define study subjects/parameters for this pipeline run
    x = xASL_init_DefineStudyData(x);

    
    %% Remove lock dirs from previous runs, if ExploreASL is not running in parallel
    if nWorkers==1
        x = xASL_init_RemoveLockDirs(x);
    end


    %% Define & print settings
    x = xASL_init_PrintCheckSettings(x);
    x = xASL_init_FileSystem(x);
end


end


%% ==================================================================================
%% ==================================================================================
function x  = xASL_init_Toolboxes(x)
%xASL_init_Toolboxes Check & load ancillary toolboxes, versions and paths.

    x.SPMDIR = fullfile(x.MyPath, 'External', 'SPMmodified');
    x.SPMpath = x.SPMDIR;
    x.SPMVERSION = 'SPM12';

    if isdeployed
        spm('asciiwelcome');
        spm('defaults','pet');
        spm_jobman('initcfg');
    else
        addpath(fullfile(x.SPMpath ,'compat'));
    end

    spm_get_defaults('cmd_line',true);

end

%% ==================================================================================
%% ==================================================================================
function [x] = xASL_init_RemoveLockDirs(x)
%xASL_init_RemoveLockDirs Remove 'lock-dir' if present from aborted previous run, for current subjects only


    % LockDir within 2 directories (e.g. T1, FLAIR or ASL)
    LockDir = fullfile(x.D.ROOT, 'lock');

    if exist(LockDir, 'dir')
%             fprintf('%s\n','Searching for locked previous ExploreASL image processing');
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
%                 fprintf('%s\n', 'No locked folders found from previous ExploreASL image processing');
        end
    end

end


%% ==================================================================================
%% ==================================================================================
function [x] = xASL_init_DefinePaths(x)
%xASL_init_DefinePaths Define paths used by ExploreASL

%% General

%% Prefixes standard space
x.D.CBFPreFix_Resliced  = 'qCBF';
x.D.c_PreFix{1} = 'rc1T1';
x.D.c_PreFix{2} = 'rc2T1';

% Atlases and templates
if isfield(x, 'MyPath')
    x.D.MapsDir             = fullfile(x.MyPath, 'Maps');
    x.D.MapsSPMmodifiedDir  = fullfile(x.MyPath, 'External', 'SPMmodified', 'MapsAdded');
    x.D.ResliceRef          = fullfile(x.MyPath, 'External', 'SPMmodified', 'MapsAdded', 'rgrey.nii');
    x.D.IdentityTransfRef   = fullfile(x.MyPath, 'External', 'SPMmodified', 'MapsAdded', 'Identity_Deformation_y_T1.nii');
    x.D.TemplateDir         = fullfile(x.MyPath, 'Maps', 'Templates');
    x.D.AtlasDir            = fullfile(x.MyPath, 'External', 'AtlasesNonCommercial');
else
    warning('MyPath field not defined...');
end

%% Study-specific
if and(isfield(x.D, 'ROOT'), isfield(x, 'ProcessData'))
    if x.ProcessData
        x.D.PopDir = fullfile(x.D.ROOT,'Population');

        % Structural module
        x.D.T1CheckDir          = fullfile(x.D.PopDir, 'T1Check');
        x.D.TissueVolumeDir     = fullfile(x.D.PopDir, 'TissueVolume');
        x.D.CoregDir            = fullfile(x.D.PopDir, 'T1wCoregCheck');
        x.D.FLAIR_CheckDir      = fullfile(x.D.PopDir, 'FLAIRCheck' );
        x.D.FLAIR_REGDIR        = fullfile(x.D.PopDir, 'FLAIRReg'   );
        x.D.FlowFieldCheck      = fullfile(x.D.PopDir, 'FlowFieldCheck' );
        x.D.LongRegCheckDir     = fullfile(x.D.PopDir, 'LongRegCheck');
        x.D.LesionCheckDir      = fullfile(x.D.PopDir, 'LesionCheck');
        x.D.ROICheckDir         = fullfile(x.D.PopDir, 'ROICheck');

        % ASL module
        x.D.ASLCheckDir         = fullfile(x.D.PopDir, 'ASLCheck');
        x.D.MotionDir           = fullfile(x.D.PopDir, 'MotionASL');
        x.D.ExclusionDir        = fullfile(x.D.PopDir, 'Exclusion');
        x.D.DICOMparameterDir   = fullfile(x.D.PopDir, 'DICOMparameters');
        x.D.SNRdir              = fullfile(x.D.PopDir, 'SD_SNR');
        x.D.M0CheckDir          = fullfile(x.D.PopDir, 'M0Check');
        x.D.M0regASLdir         = fullfile(x.D.PopDir, 'M0Reg_ASL');
        x.D.SliceCheckDir       = fullfile(x.D.PopDir, 'SliceGradientCheck');
        x.D.RawDir              = fullfile(x.D.PopDir, 'RawSourceIMCheck');
        x.D.RawEPIdir           = fullfile(x.D.PopDir, 'Raw_EPI_Check');
        x.D.T1_ASLREGDIR        = fullfile(x.D.PopDir, 'T1_ASLReg');
        x.D.TTCheckDir          = fullfile(x.D.PopDir, 'ATT_Check');
        x.D.TemplatesStudyDir   = fullfile(x.D.PopDir, 'Templates');

        % ANALYZE module
        x.SpaghettiDir          = fullfile(x.D.PopDir, 'SpaghettiPlots');
        x.S.StatsDir            = fullfile(x.D.PopDir, 'Stats');
        x.HistogramDir          = fullfile(x.D.PopDir, 'Histograms');
        x.StatsMaps             = fullfile(x.D.PopDir, 'StatsMaps');

        xASL_adm_CreateDir(x.D.PopDir);
    end
end


end

%% ==================================================================================
%% ==================================================================================
function [x] = xASL_init_LoadDataParameterFile(x, DataParPath, SelectParFile, bUseGUI)
%xASL_init_LoadDataParameterFile

    if SelectParFile
        if bUseGUI
            if ismac % somehow macOS doesnt show the title of the question dialogue
                fprintf('Select the study-specific parameter file --->>> DataParameters.(m|json)');
            end
            [name, pathstr] = uigetfile('*.m;*.mat;*.json', 'Select the study-specific parameter file --->>> DataParameters.(m|json)');
            if sum(pathstr==0) || ~(strcmp(name(end-1:end),'.m') || strcmp(name(end-3:end),'.mat') || strcmp(name(end-4:end),'.json'))
                return
            end
            DataParPath = fullfile(pathstr,name);
        else
            fprintf('ExploreASL requires a data parameter file, which can be either a .m, .mat or .json file\n');
            DataParPath = input('Please provide the full path of the DataParameter file, including ": ');
        end
    end

    [pathstr, ~, Dext] = fileparts(DataParPath);
    xBackup = x;

    if strcmp(Dext,'.mat') % Now reading a mat-file, soon to be .json table complying with BIDS
        %% Mat file can only contain single variable, that should contain the x/parameters
        TempVar = load(DataParPath);
        FieldN = fields(TempVar);
        x = TempVar.(FieldN{1});
    elseif strcmp(Dext,'.json')
        % JSON import
        try
            x = xASL_import_json(DataParPath);
        catch ME1
            % if this fails, try to recreate the json file from an .m file,
            % if it exists
            [Fpath, Ffile] = fileparts(DataParPath);
            DataParPath = fullfile(Fpath, [Ffile '.m']);
            if exist(DataParPath, 'file')
                warning('Invalid DataPar JSON file, trying to repair from detected .m file');
                fprintf('A common issue is needing escaping e.g. "\\d" instead of "\d"\n');
                try
                    PathJSON = xASL_init_ConvertM2JSON(DataParPath);
                    DataParPath = PathJSON;
                    x = xASL_import_json(DataParPath);
                catch ME2
                    fprintf('%s\n', ME1.message);
                    fprintf('%s\n', ME2.message);
                    fprintf('A common issue is needing escaping e.g. "\\d" instead of "\d"\n');
                    error('Something went wrong loading the DataParFile');
                end
            else
                warning('Invalid DataPar JSON file');
                fprintf('A common issue is needing escaping e.g. "\\d" instead of "\d"\n');
                fprintf('%s\n', ME1.message);
                error('Something went wrong loading the DataParFile');
            end
        end
    elseif strcmp(Dext,'.m')
        try
            %% Backward compatibility
            PathJSON = xASL_init_ConvertM2JSON(DataParPath); % convert .m to .json
            x = xASL_import_json(PathJSON);
        catch ME1
            try
                % Bypass eval error stuff with long names, spaces etc
                %% BACKWARD COMPATIBILITY FOR NOW, REPLACE WITH THE xASL_init_ConvertM2JSON ABOVE
                TempDataParPath = fullfile(pathstr,'TempDataPar.m');
                copyfile(DataParPath, TempDataParPath,'f' ); % overwrite ,if exist

                pathstr = fileparts(TempDataParPath);
                BackupCD = pwd;
                cd(pathstr);
                x = TempDataPar;
                cd(BackupCD);
                if exist(TempDataParPath,'file')
                    delete(TempDataParPath);
                end
            catch ME2
                fprintf('%s',ME1.message);
                fprintf('%s',ME2.message);
                error('Something went wrong loading the DataParFile');
            end
        end
    end
    % Put x fields back from backup
    FieldsAre = fields(xBackup);
    for iField=1:length(FieldsAre)
        if ~isfield(x,(FieldsAre{iField}))
            x.(FieldsAre{iField}) = xBackup.(FieldsAre{iField});
        end
    end

    if ~isfield(x,'D')
        x.D = struct;
	end

    if ~isfield(x.D,'ROOT')
		if isfield(x,'ROOT')
			x.D.ROOT = x.ROOT;
		else
			x.D.ROOT = pathstr; % default
		end
    end

    if ~exist(x.D.ROOT, 'dir')
        warning([x.D.ROOT ' didnt exist as folder, trying path of DataPar file']);
        x.D.ROOT = pathstr;
        x.ROOT = pathstr;
    end

end

%% ==================================================================================
%% ==================================================================================
function [x] = xASL_init_DefineIndependentSettings(x)
%xASL_init_DefineIndependentSettings Define ExploreASL environment
%parameters, independent of loaded data



%% --------------------------------------------------------------------------------------------------------------------
%% Default parameters
x.stopaftererrors = Inf; % set to a high number (or Inf) when running a complete batch overnight
x.dryrun = false; % set to true to skip all real processing & just print all parameters
x.bOverwrite = true;

if ~exist('groot','builtin')
    % before R2012b
    set(0,'DefaultTextInterpreter','none')
else
    % R2012b+
    set(groot,'DefaultTextInterpreter','none')
end

% Get version
if ~isdeployed
    VersionPath = xASL_adm_GetFileList(x.MyPath, '^VERSION.*$', 'FPList', [0 Inf]);
    if isempty(VersionPath)
        warning('Could not obtain ExploreASL version, version file missing');
    else
        [~, Fname, Fext] = fileparts(VersionPath{1});
        x.Version = [Fname(9:end) Fext];
    end
else
    % Output of compiled ExploreASL version
    try
        versionFile = dir(fullfile(x.MyPath, 'VERSION*'));
        versionFile = versionFile.name;
        x.Version = versionFile(9:end);
    catch
        x.Version = '1.0.0 (TMP)';
    end
end

if ~isfield(x,'Q')
    x.Q = struct;
end

x = xASL_init_DefinePaths(x);
x = xASL_init_Toolboxes(x); % Initialize toolboxes
x = xASL_init_VisualizationSettings(x); % visual settings



end


%% ==================================================================================
%% ==================================================================================
function [x] = xASL_init_DefineDataDependentSettings(x)
%xASL_init_DefineDataDependentSettings Define ExploreASL environment
%parameters, dependent of loaded data


%% --------------------------------------------------------------------------
%% Reproducibility testing
if ~isfield(x,'bReproTesting')
    x.bReproTesting = false;
end

% If the reproducibility is on, then delete the old RMS file
if isfield(x, 'bReproTesting')
    if x.bReproTesting
        xASL_delete(fullfile(x.D.ROOT, 'RMS_Reproducibility.mat'))
    end
end

%% --------------------------------------------------------------------------
%% Setting the option for pediatric template (this is normally set only for specific dataset and general xASL initialization does not have it)
% Set if the pediatric template field is set correctly
if ~isfield(x,'Pediatric_Template') || isempty(x.Pediatric_Template)
    x.Pediatric_Template = false;
end

% We need to define which of the pediatric templates to use
if x.Pediatric_Template
	% Set the default
	if ~isfield(x,'Pediatric_Type') || isempty(x.Pediatric_Type)
		x.Pediatric_Type = 'infant-1yr';
	end
	
	switch (x.Pediatric_Type)
		case {'1yr','infant_1yr','Infant_1yr'}
			x.Pediatric_Type = 'infant-1yr';
		case {'2yr','infant_2yr','Infant_2yr'}
			x.Pediatric_Type = 'infant-2yr';
		case {'neo','infant_neo','Infant_neo','neonate'}
			x.Pediatric_Type = 'infant-neo';
	end
end


%% --------------------------------------------------------------------------
%% Atlases and templates in pediatric version
if x.Pediatric_Template
	x.D.MapsDir             = fullfile(x.MyPath,'Maps', x.Pediatric_Type);
	x.D.MapsSPMmodifiedDir  = fullfile(x.MyPath,'External', 'SPMmodified', 'MapsAdded', x.Pediatric_Type);
	x.D.ResliceRef          = fullfile(x.MyPath,'External', 'SPMmodified', 'MapsAdded', x.Pediatric_Type,'rgrey.nii');
	x.D.IdentityTransfRef   = fullfile(x.MyPath,'External', 'SPMmodified', 'MapsAdded', x.Pediatric_Type,'Identity_Deformation_y_T1.nii');
	x.D.TemplateDir         = fullfile(x.MyPath,'Maps', 'Templates', x.Pediatric_Type);
	x.D.AtlasDir            = fullfile(x.MyPath,'External', 'AtlasesNonCommercial', x.Pediatric_Type);
end

if ~isfield(x, 'SegmentSPM12') && isfield(x, 'Segment_SPM12')
    warning('Please use input parameter SegmentSPM12 instead of Segment_SPM12 (legacy)');
    fprintf(['Using legacy option: x.Segment_SPM12 = ' xASL_num2str(x.Segment_SPM12) '\n']);
    x.SegmentSPM12 = x.Segment_SPM12;
end

%% --------------------------------------------------------------------------
%% Manage input parameters ExploreASL course
Fields = {'bLesionFilling' 'bAutoACPC' 'SegmentSPM12' 'M0_conventionalProcessing' 'bGetControlLabelOrder'};
Defaults = [true true false false true];

for iL=1:length(Fields)
    if ~isfield(x,Fields{iL})
        x.(Fields{iL}) = Defaults(iL);
    elseif isempty(x.(Fields{iL}))
        x.(Fields{iL}) = Defaults(iL);
    elseif ~islogical(x.(Fields{iL}))
        if x.(Fields{iL})==1
            x.(Fields{iL}) = true;
        elseif x.(Fields{iL})==0
            x.(Fields{iL}) = false;
        else
            warning([Fields{iL} ' was not true or false, set to default:' num2str(Defaults(iL))]);
            x.(Fields{iL}) = Defaults(iL);
        end
    end
end


end



%% ==================================================================================
%% ==================================================================================
function x = xASL_init_PrintCheckSettings(x)
%xASL_init_PrintCheckSettings Check whether pre-defined settings existed in DATA_PAR.m
% Prints these on the screen as the start of the pipeline.
% Runs following steps:
% 1) Set default settings if not defined
% 2) Print data/study specific settings
% 3) Print warnings

%% -----------------------------------------------------------------------
%% 1) Set default settings if not defined
if ~isfield(x,'Quality') || (x.Quality~=0 && x.Quality~=1)
    x.Quality = 1;
    fprintf('%s\n', 'Default Quality=1 used (optimal quality)');
end
if ~isfield(x,'DELETETEMP') || (x.DELETETEMP~=0 && x.DELETETEMP~=1)
    x.DELETETEMP = 1;
%     fprintf('%s\n','Default x.DELETETEMP=1 used (delete files temporarily used for processing)');
end

%% -----------------------------------------------------------------------
%% 2) Print data/study specific settings
fprintf('%s\n','-------------------------------------------');
fprintf('%s\n\n','ExploreASL will run with following settings:');
fprintf('%s\n\n',['Root folder = ' x.D.ROOT]);

if x.nWorkers>1
    fprintf(['I am worker ' num2str(x.iWorker) '/' num2str(x.nWorkers) '\n']);
    fprintf('Note that the resulting number of scans mentioned below applies only to this worker\n');
end

fprintf('%s\n',[num2str(x.nTotalSubjects) ' scans - ' num2str(x.nExcluded) ' exclusions, resulting in ' num2str(x.nSubjects) ' scans of: ']);

for iT=1:x.nTimePointsTotal
    fprintf('%s\n',['Longitudinal timePoint ' num2str(iT) ' = ' num2str(x.nTimePointTotalSubjects(iT)) ' scans - ' num2str(x.nTimePointExcluded(iT)) ' exclusions = ' num2str(x.nTimePointSubjects(iT)) ' scans']);
end

fprintf('%s\n',['ASL sessions: ' num2str(x.nSessions)]);

fprintf('\n%s\n','Ancillary data, sets:');
if isfield(x.S,'SetsID')
        fprintf('%s\n',[num2str(size(x.S.SetsID,2)) ' sets are defined for ' num2str(size(x.S.SetsID,1)) ' "SubjectsSessions":']);

        for iSet=1:size(x.S.SetsID,2)
            fprintf(['Set ' num2str(iSet) ' = "' x.S.SetsName{iSet} '" options ']);
            for iOption=1:size(x.S.SetsOptions{iSet},2)
                fprintf(['"' x.S.SetsOptions{iSet}{iOption} '"']);
                if iOption~= size(x.S.SetsOptions{iSet},2)
                    fprintf(' & ');
                end
            end
            if      x.S.Sets1_2Sample(iSet)==1
                    fprintf(', codes for paired data');
            elseif  x.S.Sets1_2Sample(iSet)==2
                    fprintf(', codes for two-sample data');
            elseif  x.S.Sets1_2Sample(iSet)==3
                    fprintf([', continuous variate (with ' num2str(length(unique(x.S.SetsID(:,iSet)))) ' unique values)']);
            end                    

            fprintf('\n');
        end
else    
    fprintf('%s\n','No sets are defined');
end

if ~isfield(x,'M0')
%     warning('M0 option missing!');
else
    fprintf('\n%s\n',['M0 option selected is "' num2str(x.M0) '"']);
end

fprintf('%s\n',['x.DELETETEMP = ' num2str(x.DELETETEMP) ' (delete temporary files)']);
fprintf('%s\n',['x.Quality    = ' num2str(x.Quality) ' (0 = fast try-out; 1 = normal high quality)']);


%% -----------------------------------------------------------------------
%% 3) Print warnings
fprintf('\n%s\n\n','---------------------------------------------');

field_symbol = {'subject_regexp'};

for iField=1:length(field_symbol)
    if ~isfield(x,field_symbol{iField})
        warning(['x.' field_symbol{iField} ' was not defined in DATA_PAR.m!'])
    end
end

if ~isfield(x,'D')
    warning('x.D didn''nt exist');
else
    field_symbol = {'ROOT'};
    for iField=1:length(field_symbol)
        if ~isfield(x.D,field_symbol{iField})
            warning(['x.D.' field_symbol{iField} ' was not defined in DATA_PAR.m!'])
        end
    end
end

if ~isempty(regexp(x.subject_regexp, '^(\^|)\.\*(\$|)$'))
    warning('Subject regexp not specific! Check that no wrong folders are included as subjects');
end

fprintf('\n');

end

