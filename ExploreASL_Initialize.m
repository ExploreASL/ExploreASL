function [x] = ExploreASL_Initialize(DataParPath, ProcessData, iWorker, nWorkers)
%ExploreASL_Initialize Initializes ExploreASL
%
% FORMAT: [x] = ExploreASL_Initialize([DataParPath, ProcessData, iWorker, nWorkers])
%
% INPUT:
%   DataParPath - path to data parameter file (OPTIONAL, required when ProcessData=true, will then be prompted if omitted)
%   ProcessData - true if running pipeline on data, false if only initializing ExploreASL (e.g. path management etc, REQUIRED, will be prompted if omitted)
%   iWorker     - allows parallelization when called externally. iWorker defines which of the parallel ExploreASL calls we are (OPTIONAL, DEFAULT=1)
%   nWorkers    - allows parallelization when called externally. nWorkers defines how many ExploreASL calls are made in parallel (OPTIONAL, DEFAULT=1)
%
% OUTPUT:
%   x           - struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This initialization wrapper initializes ExploreASL:
% managing paths, deployment,
%
% xASL_init_ExploreASL_directories - manages folders for ExploreASL, and sets and creates study/data-folders
% xASL_init_Toolboxes              - initialization third-party toolboxes, e.g. SPM, dip_image (soon to be removed)
% xASL_init_VisualizationSettings  - defines visualization settings for
%                                    visual QC figure printing (type help xASL_init_VisualizationSettings for more information)
% xASL_init_PrintCheckSettings     - prints summarized data parameters and warnings
% xASL_init_FileSystem             - dirty initialization of common filenames used throughout the pipeline
% xASL_init_PopulationSettings     - defines additional visualization settings (to be merged with xASL_init_VisualizationSettings)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE for GUI: ExploreASL_Initialize
% EXAMPLE for calling externally: ExploreASL_Initialize('//MyDisk/MyStudy/DataPar.m', true, 1, 1);
% EXAMPLE for debugging/initialization only: [x] = ExploreASL_Initialize('',0);
% __________________________________
% Copyright 2015-2019 ExploreASL

%% Load with or without dataset
% Construct a questdlg with three options

% Defaults
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
    if ProcessData==1 || ProcessData==0
        x.ProcessData = ProcessData;
    else
        warning('argin "ProcessData" should be 1 or 0 (true or false)');
    end
end

%% Find out whether we are run in GUI or CLI
 if ~usejava('desktop') || ~usejava('jvm') || ~feature('ShowFigureWindows')
     UseGUI = false;
 else
     UseGUI = true;
 end

if ~isfield(x,'ProcessData')

    x.InitChoice = 'Process dataset';
    if ~exist('DataParPath','var')
        ParFileExists = false;
    elseif isempty(DataParPath)
        ParFileExists = false;
    else
        ParFileExists = true;
    end

    if ~ParFileExists
        if UseGUI
            x.InitChoice = questdlg('Would you like to load a dataset or only initialize ExploreASL (set paths etc)?', ...
                'Initialize ExploreASL', 'Process dataset','Only initialize ExploreASL','Cancel','Cancel');
        else
            fprintf('Would you like to load a dataset or only initialize ExploreASL (set paths etc)?\n');
            cliChoice = input('Please press [1] for "Process dataset", [2] for "Only initialize ExploreASL", [3] to cancel: ');

            switch cliChoice
                case 1
                    x.InitChoice = 'Process dataset';
                case 2
                    x.InitChoice = 'Only initialize ExploreASL';
                case 3
                    x.InitChoice = 'Cancel';
                otherwise
                    fprintf('Wrong choice, please choose 1 2 or 3: ');
                    return;
            end
        end
    end
    % Handle response
    switch x.InitChoice
        case 'Process dataset'
            fprintf('%s\n','Processing & analyzing dataset')
            x.ProcessData = true;
        case 'Only initialize ExploreASL'
            fprintf('%s\n','Only initializing ExploreASL (set paths etc)');
            x.ProcessData = false;
        case 'Cancel'
            fprintf('%s\n','Exiting ExploreASL, please ignore the errors');
            x.ProcessData = false;
            error('Canceled ExploreASL');
        otherwise
            x.ProcessData = false;
            error('Something wrong, please retry');
    end
end

%% -----------------------------------------------------------------------------
%% Get ExploreASL path

% Check if the current directory is the ExploreASL directory
CurrCD = pwd;
if exist(fullfile(CurrCD, 'ExploreASL_Master.m'), 'file')
    MyPath2  = CurrCD;
end

% Check whether MyPath is correct, otherwise obtain correct folder
if ~isfield(x, 'MyPath')
    x.MyPath  = '/DummyPath';
end

if exist('MyPath2','var')
    x.MyPath = MyPath2;
end

MasterScriptPath = fullfile(x.MyPath, 'ExploreASL_Master.m');

if ~isdeployed
    if ~exist(MasterScriptPath,'file')
        if UseGUI
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
    disp('------------ ctfroot -------------');
    disp(ctfroot);
    disp('------------ x.MyPath ------------');
    disp(x.MyPath);
    disp('----------------------------------');
    % pause;
end

% Go to ExploreASL folder
cd(x.MyPath);


%% Add paths
if ~isdeployed
% Define paths (should be equal when loading data or only initializing
% ExploreASL
    addpath(x.MyPath);

	%fullfile('Development','dicomtools'), ...
    subfolders_to_add = {   'Functions', ...
                            'mex', ...
                            'Modules', ...
                            fullfile('Modules', 'SubModule_Structural'),...
                            fullfile('Modules', 'SubModule_ASL'),...
                            fullfile('Modules', 'SubModule_Population'),...
                            'Development' ...
                            'External' ...
                            fullfile('External','DCMTK'), ...
                            fullfile('External','ExploreQC'), ...
							fullfile('External','SPMmodified'), ...
                            fullfile('External','SPMmodified','xASL'),...                            
                            fullfile('External','SPMmodified','toolbox','cat12'), ...
							fullfile('External','SPMmodified','toolbox','LST'), ...
							fullfile('External','SPMmodified','toolbox','OldNorm')};

    for ii=1:length(subfolders_to_add)
        addpath(fullfile(x.MyPath,subfolders_to_add{ii}));
    end
end

if x.ProcessData
    %% Load dataset parameter file

    if SelectParFile
        if UseGUI
            [name, pathstr] = uigetfile('*.m;*.mat;*.json', 'Select the study-specific PARAMETER file --->>> DATA_PAR.m');
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
        x = xASL_import_json(DataParPath);
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
    x.ProcessData = xBackup.ProcessData;
    x.iWorker = xBackup.iWorker;
    x.nWorkers = xBackup.nWorkers;

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





%% -----------------------------------------------------------------------------
%% Initialization

x.stopaftererrors     = Inf;      % set to a high number (or Inf) when running a complete batch overnight
x.dryrun              = false;    % set to true to skip all real processing & just print all parameters

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
    x.Version = ' compiled version, check in filename';
end

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
%xASL_adm_CheckPermissions(x.MyPath, false); % doesn't have to be executed


%% -----------------------------------------------------------------------------
%% Common settings and definitions
x                     = xASL_init_ExploreASL_directories(x);
x.OVERWRITE           = true;

x                     = xASL_init_Toolboxes(x);
[x]                   = xASL_init_VisualizationSettings(x); % visual settings


%% Reproducibility testing
if ~isfield(x,'bReproTesting')
    x.bReproTesting = false;
end

if ~x.ProcessData
    fprintf('%s\n','ExploreASL initialized <--');
else

    if isempty(x.D.ROOT)
        error('No dir defined for data');
    end

    % Fix a relative path
    if strcmp(x.D.ROOT(1), '.')
        cd(x.D.ROOT);
        x.D.ROOT = pwd;
    end

    %% Define subjects/parameters
    x = xASL_init_DefineSets(x);

    
    
    
    
    
    
    
    
    if ~isfield(x,'name'); x.name = ''; end


    %% Remove 'lock-dir' if present from aborted previous run, for current subjects only
    % LockDir within 2 directories (e.g. T1, FLAIR or ASL)
    LOCKDIR = fullfile(x.D.ROOT,'lock');
    % COMMENTED THIS OUT FOR PARALLELIZATION
%     if  isdir(LOCKDIR)
%         fprintf('%s\n','Searching for locked previous ExploreASL image processing');
%         LockDirFound    = 0;
%         [ LockDir, ~] = xASL_adm_FindByRegExp( fullfile(x.D.ROOT,'lock'), {'(ASL|Struct|LongReg_T1)',x.subject_regexp,'.*module.*','^(locked)$'}, 'Match','Directories');
%         if  ~isempty(LockDir)
%             for iL=1:length(LockDir)
%                 if isdir(LockDir{1})
%                     fprintf('%s\n',[LockDir{1} ' detected, removing']);
%                     rmdir(LockDir{1},'s');
%                 end
%             end
%             LockDirFound    = 1;
%         end
%
%         % LockDir within 2 directories (e.g. DARTEL)
%         [ LockDir, ~] = xASL_adm_FindByRegExp( fullfile(x.D.ROOT,'lock'), {'(Population|DARTEL_T1)','.*module.*','^(locked)$'}, 'Match','Directories');
%         if  ~isempty(LockDir)
%             for iL=1:length(LockDir)
%                 fprintf('%s\n',[LockDir{1} ' detected, removing']);
%                 rmdir(LockDir{1},'s');
%             end
%             LockDirFound    = 1;
%         end
%
%         if  LockDirFound==0
%             fprintf('%s\n','No locked dirs found from previous ExploreASL image processing');
%         end
%     end

    %% Print settings to check
    x = xASL_init_PrintCheckSettings( x );
    x = xASL_init_FileSystem(x);
    x = xASL_init_PopulationSettings( x);
end

%% Manage input parameters ExploreASL course
Fields   = {'bLesionFilling' 'bAutoACPC' 'Segment_SPM12' 'bRegistrationContrast' 'M0_conventionalProcessing' 'bGetControlLabelOrder'};
Defaults = [true true false true false true];

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

% If the reproducibility is on, then delete the old RMS file
if isfield(x,'bReproTesting')
    if x.bReproTesting
        xASL_delete(fullfile(x.D.ROOT,'RMS_Reproducibility.mat'))
    end
end

end
