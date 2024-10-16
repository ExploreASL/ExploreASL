function [comparison] = xASL_test_TestExploreASL(TestDirOrig, TestDirDest, RunMethod, bTestSPM, MatlabPath, EmailAddress, Password, bOverwrite, testDataUsed, RunTimePath, bPull)
%xASL_test_TestExploreASL Do a thorough test of the validity and reproducibility of ExploreASL
%
% FORMAT: [ResultsTable] = xASL_test_TestExploreASL(TestDirOrig, TestDirDest, RunMethod, bTestSPM, 
%                          MatlabPath, EmailAddress, Password, bOverwrite, testDataUsed, RunTimePath, bPull)
% 
% INPUT:
%   TestDirOrig - path to root folder containing all datasets to test processing on (REQUIRED)
%   TestDirDest - path to folder where the results are stored (this is temporarily, automatically deleted upon restart, REQUIRED)
%   RunMethod   - Choose a value for how to test ExploreASL (REQUIRED)
%                 Option 1 = run ExploreASL serially
%                 Option 2 = run ExploreASl parallel (start new MATLAB instances)
%                 FUTURE Option 3 = run ExploreASL compilation serially
%                 FUTURE Option 4 = run ExploreASL compilation parallel
%   bTestSPM    - boolean for testing if SPM standalone with xASL modifications works (DEFAULT=true)
%   MatlabPath  - path to matlab executable or compilation bash script (OPTIONAL, required in some cases)
%   EmailAddress- string with e-mail address for gmail account to use (OPTIONAL, DEFAULT = skip e-mailing results)
%   Password    - string with password for this gmail account (REQUIRED when EmailAddress provided)
%   bOverwrite  - Overwrite existing test results (OPTIONAL, DEFAULT=true);
%   testDataUsed- Option 1: ExploreASL/External/TestDataSet as an input
%                 Option 0: Other directory used as TestDirOrig (DEFAULT)
%   RunTimePath - When using a compiled version, the location of the
%                 Matlab RunTime libraries (e.g. '/usr/local/MATLAB/MATLAB_Runtime/v96')
%   bPull       - pull new version of the software (OPTIONAL, DEFAULT=true)
%                 
% OUTPUT:
%   ResultsTable - Table containing all results from the test runs
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function will run ExploreASL on several different
%              datasets, do perform a full and thorough test of ExploreASL.
%              Different environment variables include:
%              x.settings.Quality 0 and 1
%              
%              Different data setups include:
%              ASL readouts (3D spiral, 3D GRASE, 2D EPI)
%              ASL Manufacturers (GE, Philips, Siemens)
%              With/without background suppression
%              With/without FLAIR processing (LST LGA or LPA)
%              With/without lesion masking from tumor
%              Performed as first or second run of ExploreASL, fully or partly done before
%              Longitudinal or cross-sectional data
%              FEAST
%              With/without disabling quantification
%              With/without saving 4D PWI/CBF maps
%              With/without M0
%              With single or multiple ASL sessions (i.e. runs)
%
%              This function performs the following steps:
%
%              1. Pull latest GitHub version
%              2. Copy all data for testing
%              3. Initialize & test standalone SPM on low quality
%              4. Test ExploreASL itself
%              5. Pause until all results exist (if running parallel in background)
%              6. Compile results table
%              7. Compare table with reference table
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:
%
% Jan:           [ResultsTable] = xASL_test_TestExploreASL('/pet/projekte/asl/data/ExploreASL_TestCases', '/pet/projekte/asl/data/ExploreASL_TempRes', 1);
% Henk on MacOS: [ResultsTable] = xASL_test_TestExploreASL('/Users/henk/surfdrive/HolidayPics/ExploreASL_TestCases', '/Users/henk/ExploreASL/ASL/ExploreASL_TestCasesProcessed', 1, 0,[],'henkjanmutsaerts@gmail.com');
% VUmc server:   [ResultsTable] = xASL_test_TestExploreASL('/radshare/ExploreASL_Test/ExploreASL_TestCases', '/radshare/ExploreASL_Test/ExploreASL_TestCasesProcessed', 1);
% __________________________________
% Copyright (c) 2015-2022 ExploreASL

% ============================================================
%% Admin

% Run ExploreASL to get directories
if isempty(which('ExploreASL'))
    cd ..;
else
    cd(fileparts(which('ExploreASL')));
end

% Validate the input options
if nargin < 3 || isempty(TestDirOrig) || isempty(TestDirDest) || isempty(RunMethod)
	error('Require three input parameters -TestDirOrig, TestDirDest, RunMethod');
end

if nargin < 4 || isempty(bTestSPM)
	bTestSPM = true;
end

if nargin < 5 || isempty(MatlabPath)
	MatlabPath = 'matlab';
end

if nargin < 6 || isempty(EmailAddress)
	EmailAddress = '';
end

if nargin < 7 || isempty(Password)
	Password = '';
end

if nargin < 8 || isempty(bOverwrite)
	bOverwrite = true;
end

if nargin < 9 || isempty(testDataUsed)
	testDataUsed = 0;
end

if nargin < 10 || isempty(RunTimePath)
	RunTimePath = '';
end

if nargin < 11 || isempty(bPull)
	bPull = true;
end

% Check Matlab path and Runtime path for corresponding run method
if RunMethod>2
    if isempty(MatlabPath) || ~exist(MatlabPath, 'file') || ~strcmp(MatlabPath(end-2:end),'.sh')
        warning('Please provide the path to the bash script calling the compiled ExploreASL, skipping');
        return;
    elseif isempty(RunTimePath) || ~exist(RunTimePath, 'dir')
        warning('Please provide the path to the Matlab Runtime installation, skipping');
        return;        
    end
end

% ============================================================
%% 1) Pull latest GitHub version
xASL_adm_BreakString('1. Update ExploreASL','=');

% Assuming we are in ExploreASL folder
if bPull
    xASL_system('git fetch');
    xASL_system('git pull');
end

% Initialize ExploreASL
x = ExploreASL;

% ============================================================
%% 2) Copy all data for testing
xASL_adm_BreakString('2. Copy the test data','=');

% Ask for directories if they were not defined
xASL_test_CopyTestData(TestDirDest, TestDirOrig, bOverwrite);

% ============================================================
%% 3) Initialize & test standalone SPM on low quality
xASL_adm_BreakString('3. Initialize & test standalone SPM','=');
if bTestSPM
    xASL_test_SPM(TestDirDest, testDataUsed);
end

% ============================================================
%% 4) Test ExploreASL itself

% Here we return the ExploreASL paths, which we removed above for testing SPM
x = ExploreASL;
Dlist = xASL_adm_GetFileList(TestDirDest,'^.*$','List',[0 Inf], true);
xASL_adm_BreakString('4. Test ExploreASL','=');
LogFiles = xASL_test_TestAllTestdatasets(TestDirDest, RunMethod, MatlabPath, RunTimePath, x, Dlist);

% ============================================================
%% 5) Pause until all results exist (if running parallel in background)
xASL_adm_BreakString('5. Pause','=');
if RunMethod==2 || RunMethod==4
    CountTime = 0;
    TimeStepSeconds = 30;
    fprintf(xASL_adm_ConvertSeconds2TimeString(CountTime));

    while max(cellfun(@(y) ~exist(y,'file'), LogFiles)) % logs don't exist
        pause(TimeStepSeconds);
        CountTime = CountTime+TimeStepSeconds;
        TimeString = xASL_adm_ConvertSeconds2TimeString(CountTime);
        fprintf(['\b\b\b\b\b\b' TimeString]);
    end
    fprintf('\n');
end

% ============================================================
%% 6) Compile results table
xASL_adm_BreakString('6. Compile results table and compare with references','=');
comparison = xASL_test_CompareReference(fullfile(x.opts.MyPath,'Testing','Reference','ReferenceValues.tsv'), TestDirDest, TestDirOrig);

% Comparison with mat file
try
    % Determine the difference table
    comparison.DifferenceTable = xASL_test_DetermineDifferenceTable(TestDirOrig, comparison.ResultsTable, comparison.ResultTableFile);
    
    % E-Mail the results
    if ~isempty(EmailAddress)
        xASL_test_EmailResults(EmailAddress, Password, comparison.DifferenceTable);
    end
catch ME
    warning('Something went wrong in trying to create difference table & mailing it to receivers');
    fprintf('%s\n', ME.message);
end
    
save(comparison.SaveFile, 'comparison');
end

%% Copy the test data
function xASL_test_CopyTestData(TestDirDest, TestDirOrig, bOverwrite)
    if isempty(TestDirDest)
        TestDirDest = uigetdir(pwd, 'Select testing directory...');
    end
    if isempty(TestDirOrig)
        TestDirOrig = uigetdir(pwd, 'Select datasets for testing...');
    end

    % Clone testdataset repository if not detected
    if ~xASL_exist(TestDirOrig, 'dir')
        TestDirRoot = fileparts(TestDirOrig);
        TestDataSetRepository = 'https://github.com/ExploreASL/TestDataSets.git';
        fprintf('%s\n', ['TestDataSet repository not found in: ' TestDirRoot]);
        fprintf('%s\n', ['Attempting to clone: ' TestDataSetRepository]);
        xASL_system(['cd ' TestDirRoot]);
        xASL_system(['git clone ' TestDataSetRepository]);
    end

    % Remove previous results
    if bOverwrite && exist(TestDirDest,'dir')
        fprintf('Deleting previous results...\n');
        if ispc
            system(['rmdir /s /q ' TestDirDest]);
        else
            system(['rm -rf ' xASL_adm_UnixPath(TestDirDest)]);
        end
    end
    % Copy data sets into testing directory
    xASL_Copy(TestDirOrig, TestDirDest);    
end


%% Determine the difference table
function DifferenceTable = xASL_test_DetermineDifferenceTable(TestDirOrig, ResultsTable, ResultTableFile)
    % Find all result tables in directory
    AllResultTables = xASL_adm_GetFsList(fullfile(TestDirOrig),'^.+\_ResultsTable.mat$',false)';
    IndexCurrentTable = find(~isempty(strfind(AllResultTables, ResultTableFile)));
    if ~isempty(IndexCurrentTable)
        AllResultTables(IndexCurrentTable) = [];
    end
    % Previous Table name
    if ~isempty(AllResultTables)
        PreviousTableName = AllResultTables{numel(AllResultTables)};
        PreviousSaveFile = fullfile(TestDirOrig, PreviousTableName);
        PreviousTable = load(PreviousSaveFile, '-mat');
        
        clear DifferenceTable
        DifferenceTable(1:size(ResultsTable,1),1:size(ResultsTable,2)) = {''};
        DifferenceTable(1,:) = ResultsTable(1,:);
        DifferenceTable(:,1) = ResultsTable(:,1);
        for iX=2:size(ResultsTable,1)
            for iY=2:size(ResultsTable,2)
                A = xASL_str2num(ResultsTable{iX,iY});
                B = xASL_str2num(PreviousTable.ResultsTable{iX,iY});
                AsymmIndex = (A-B)/(A+B);
                DifferenceTable{iX,iY} = [xASL_num2str(AsymmIndex) '%'];
            end
        end
    else
        % Could not find previous tables
        DifferenceTable = [];
    end
end


%% E-Mail the results
function xASL_test_EmailResults(EmailAddress, Password, DifferenceTable)
    % First convert table to string to send by e-mail
    NewTable{1,1} = 'mean_qCBF_TotalGM     median_qCBF_TotalGM     median_qCBF_DeepWM     CoV_qCBF_TotalGM             GMvol                 WMvol                 CSFvol             PipelineCompleted     TC_ASL_Registration    TC_M0_Registration';
    SingleEmptyString1 = repmat(' ',[1 44]);
    SingleEmptyString2 = repmat(' ',[1 27]);
    for iX=2:size(DifferenceTable,1)
        NewTable{iX,1} = '';
        for iY=2:size(DifferenceTable,2)
            if iY<5
                NewCell = SingleEmptyString1;
            else
                NewCell = SingleEmptyString2;
            end
            
            NewCell(1:length(DifferenceTable{iX,iY})) = DifferenceTable{iX,iY};
            NewTable{iX,1} = [NewTable{iX,1} NewCell];
        end
        NewTable{iX,1} = [NewTable{iX,1} DifferenceTable{iX,1}];
    end


    % See here: https://nl.mathworks.com/help/matlab/import_export/sending-email.html
    fprintf('Sending e-mail with results\n');
    setpref('Internet', 'SMTP_Server', 'smtp.gmail.com');
    setpref('Internet', 'E_mail', EmailAddress);
    setpref('Internet', 'SMTP_Username', EmailAddress);
    setpref('Internet', 'SMTP_Password', Password);
    props = java.lang.System.getProperties;
    props.setProperty('mail.smtp.auth','true');
    props.setProperty('mail.smtp.starttls.enable', 'true');
    props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
    props.setProperty('mail.smtp.socketFactory.port','465'); % or port 587
    EmailAddresses = {'Patricia.Clement@ugent.be', 'Pieter.Vandemaele@UZGENT.be', 'j.petr@hzdr.de', 'henkjanmutsaerts@gmail.com'};
    sendmail(EmailAddresses, 'ExploreASL TestRun: %AsymmetryIndexWithTemplateResults (should be <0.01%)', NewTable);
end

%% Test all individual test datasets
function LogFiles = xASL_test_TestAllTestdatasets(TestDirDest, RunMethod, MatlabPath, RunTimePath, x, Dlist)

    % Remove lock folders, useful for rerun when debugging
    LockFolders = xASL_adm_GetFileList(TestDirDest, '(?i)^locked$', 'FPListRec', [0 Inf], true);
    if ~isempty(LockFolders)
        cellfun(@(y) xASL_delete(y), LockFolders, 'UniformOutput', false);
    end

    % Evaluate parallelization in linux
    if isunix && (RunMethod==2 || RunMethod==4)
        [Result1, Result2] = system('screen -dmS TryOut exit');
        if Result1~=0
            warning('Please install screen for testing ExploreASL parallel in a mac/linux');
            fprintf('%s\n', Result2);
            error('Skipping...');
        end
    end

    % Get list of data to test
    LogFiles = cellfun(@(y) fullfile(TestDirDest,y,'log','xASL_module_Population.log'), Dlist, 'UniformOutput',false);

    % Iterate over test datasets
    for iList=1:length(Dlist)
        
        % Get files and directories
        dirBIDS = fullfile(TestDirDest,Dlist{iList});
        
        % Useful for rerun when debugging
        xASL_delete(LogFiles{iList});
        
        % Determine shown screen name
        ScreenName = ['TestxASL_' num2str(iList)];
        
        % Check if dataset exists
        if ~isempty(dirBIDS)
            try
                xASL_test_IndividualTestdataset(RunMethod, MatlabPath, RunTimePath, x, dirBIDS, ScreenName);
            catch ME
                warning(ME.identifier, 'Something went wrong: %s', ME.message);
            end
        end
    end
end


%% Test one individual test dataset
function xASL_test_IndividualTestdataset(RunMethod, MatlabPath, RunTimePath, x, dirBIDS, ScreenName)

    % Run ExploreASL
    cd(x.opts.MyPath);

    % Prepare compilation testing
    if RunMethod>2
        [Fpath, Ffile, Fext] = fileparts(MatlabPath);
        if isunix
            CompilationString = ['cd ' Fpath ';bash ' Ffile Fext ' ' RunTimePath ' ' dirBIDS];
        else
            CompilationString = ['cd ' Fpath '; ' Ffile Fext ' ' RunTimePath ' ' dirBIDS];
        end
    end

    switch RunMethod
        case 1
            % Run ExploreASL serially (can we run screen from here? or run matlab in background, linux easy)
            ExploreASL(dirBIDS, 0,  1, false);
        case 2
            % Run ExploreASl parallel (start new MATLAB instances)
            if isunix
                ScreenString = ['screen -dmS ' ScreenName ' nice -n 10 ' MatlabPath ' -nodesktop -nosplash -r '];
                RunExploreASLString = ['"cd(''' x.opts.MyPath ''');ExploreASL(''' dirBIDS ''',0,1,0);system([''screen -SX ' ScreenName ' kill'']);"'];
            else
                ScreenString = [MatlabPath ' -nodesktop -nosplash -r '];
                RunExploreASLString = ['"cd(''' x.opts.MyPath ''');ExploreASL(''' dirBIDS ''',0,1,0);system([''exit'']);"'];
            end
            system([ScreenString RunExploreASLString ' &']);
        case 3
            % Run ExploreASL compilation serially
            system(CompilationString);
        case 4 
            % Run ExploreASL compilation parallel
            if isunix
                ScreenString = ['screen -dmS ' ScreenName ' nice -n 10 '];
            else
                ScreenString = [];
            end
            
            system([ScreenString CompilationString ' &']);
        otherwise
            fprintf(2,'Unknown run method...\n');
    end
end
