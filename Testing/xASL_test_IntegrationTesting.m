function [TestResults, CheckResults] = xASL_test_IntegrationTesting
%xASL_test_IntegrationTesting Main integration testing script
%
% INPUT:        n/a
%
% OUTPUT:       TestResults  - Test results struct
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Integration testing script. Run the full pipeline on the flavor datasets using ExploreASL_Master.
%
% EXAMPLE:      [TestResults, CheckResults] = xASL_test_IntegrationTesting;
%
% Your testConfig.json could look like this e.g.:
%
% {
%    "pathExploreASL": "...//ExploreASL",
%    "pathTest": "...//FlavorTests",
%    "cmdCloneFlavors": "git clone git@github.com:ExploreASL/FlavorDatabase.git"
% }
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


    %% Initialize Paths

    % Add testing directory
    TestingPath = mfilename('fullpath');    % Path of the current script
    filepath = fileparts(TestingPath);      % Testing directory
    xASL_dir = fileparts(filepath);         % All subfolders of xASL directory
    warning('off','all')                    % Turn of warnings (only add path related)
    addpath(genpath(xASL_dir))              % Add all scripts to path (including testing directory)
    warning('on','all')
    
    % Logging table
    loggingTable = array2table(zeros(0,2), 'VariableNames',{'message','stack'});

    %% Get paths
    if exist(fullfile(filepath,'testConfig.json'),'file')
        testConfig = spm_jsonread(fullfile(filepath,'testConfig.json'));
    else
        fprintf('Please add a testConfig.json to the Testing directory of ExploreASL...\n');
        return
    end

    % Check testConfig fields
    if ~isfield(testConfig, 'pathExploreASL')
        error('Add the pathExploreASL field to your testConfig.json...');
    end
    if ~isfield(testConfig, 'pathTest')
        error('Add the pathTest field to your testConfig.json...');
    end
    if ~isfield(testConfig, 'cmdCloneFlavors')
        error('Add the cmdCloneFlavors field to your testConfig.json...');
    end

    %% Set-up
    cd(testConfig.pathExploreASL);

    % Initialize ExploreASL
    [x] = ExploreASL_Initialize;
    
    % Clean Up
    clc

    % Clone flavors if they do not exist already
    if ~exist(fullfile(testConfig.pathTest,'FlavorDatabase'), 'dir')
        xASL_adm_CreateDir(testConfig.pathTest);
        system(testConfig.cmdCloneFlavors);
    end

    % Get the up-to-date flavor list
    flavorList = xASL_adm_GetFsList(fullfile(testConfig.pathTest,'FlavorDatabase'), '^.+$', true);
    
    % Delete existing results
    if exist(fullfile(testConfig.pathTest,'FlavorDatabaseTest'), 'dir')
        xASL_delete(fullfile(testConfig.pathTest,'FlavorDatabaseTest'), 1);
    end

    %% Iterate over flavors
    for iFlavor = 1:numel(flavorList)
        % Run individual flavor
        TestResults = xASL_test_thisFlavor(testConfig, flavorList{iFlavor});
        % Check logging
        if isfield(TestResults,'logging')
            loggingTable = xASL_test_AddLoggingEntryToTable(loggingTable,TestResults.logging);
        end
    end
    
    %% Write logging table
    xASL_test_writeLoggingTable(fullfile(testConfig.pathTest,'FlavorDatabaseTest','log.xlsx'), loggingTable);
    
    %% Compare rawdata
    for iFlavor = 1:numel(flavorList)
        % Compare rawdata of individual flavor
        CheckResults(iFlavor) = xASL_test_thisFlavorCheckRawdata(testConfig.pathTest, ...
                                                    fullfile(testConfig.pathTest,'FlavorDatabase'), ...
                                                    flavorList{iFlavor});
    end

    % Fallback
    if ~exist('TestResults','var')
        TestResults = struct;
    end
    % Fallback
    if ~exist('CheckResults','var')
        CheckResults = struct;
    end
    
    % Go back to ExploreASL directory and reset everything
    cd(testConfig.pathExploreASL);
    close all
    diary off

end


%% Run individual flavor dataset
function xFlavor = xASL_test_thisFlavor(testConfig, flavorName)

    % Get paths
    testingRoot = testConfig.pathTest;
    databaseRoot = fullfile(testConfig.pathTest,'FlavorDatabase');
    
    % Copy flavor
    fprintf('Copy %s...\n', flavorName);
    xASL_Copy(fullfile(databaseRoot,flavorName),fullfile(testingRoot,'FlavorDatabaseTest',flavorName),1);
    
    % Remove rawdata folder from testing dataset
    xASL_delete(fullfile(testingRoot,'FlavorDatabaseTest',flavorName,'rawdata'),1);
    
    % Create default dataPar.json
    defaultDataPar.x.subject_regexp = '^sub-.*$';
    defaultDataPar.x.bUseMNIasDummyStructural = 1; % when no structural data, use ASL-MNI registration
    defaultDataPar.x.Quality = 0; % speed up testing
    defaultDataPar.x.DELETETEMP = 1;
    
    % Save JSON
    spm_jsonwrite(fullfile(testingRoot,'FlavorDatabaseTest',flavorName,'dataPar.json'),defaultDataPar);
    
    % Remove Matlab paths (independent initialization)
    xASL_adm_RemoveDirectories(testConfig.pathExploreASL);

    % Run test
    try
        [xFlavor] = ExploreASL_Master(fullfile(testingRoot,'FlavorDatabaseTest',flavorName),1,[1 1 0]);
    catch ME
        fprintf(2,'%s\n', flavorName);
        fprintf(2,'Error: %s\n', ME.message);
        if size(ME.stack,1)>=1
            fprintf(2,'File: %s\n', ME.stack(1).name);
            fprintf(2,'Line: %d\n', ME.stack(1).line);
        end
    end
    
    % Fallback
    if ~exist('xFlavor','var')
        xFlavor = struct;
    end
    

end


%% Check rawdata of individual flavor dataset
function checkResults = xASL_test_thisFlavorCheckRawdata(testingRoot, databaseRoot, flavorName)

    % Compare rawdata
    fprintf('Compare rawdata of %s...\n', flavorName);
    
    % Get paths
    pathDataA = fullfile(databaseRoot,flavorName,'rawdata');
    pathDataB = fullfile(testingRoot,'FlavorDatabaseTest',flavorName,'rawdata');
    
    % Copy based comparison
    if exist(pathDataA,'dir') && exist(pathDataB,'dir')
        [identical,results] = xASL_bids_CompareStructures(pathDataA,pathDataB,1,[],[],1);
    else
        identical = 0;
        results = struct;
    end
    
    % Return results
    checkResults.name = flavorName;
    checkResults.results = results;
    checkResults.identical = identical;

end


%% Add entry to log table
function logTable = xASL_test_AddLoggingEntryToTable(logTable,logStruct)

    % Get number of log entries
    numLogEntries = size(logStruct,2);
    
    % Iterate over log entries
    for iEntry=1:numLogEntries
        thisStruct = logStruct(iEntry);
        thisStruct.stack = xASL_test_StackToString(thisStruct.stack);
        thisStruct.message = cellstr(thisStruct.message);
        thisStruct.stack = cellstr(thisStruct.stack);
        thisRow = struct2table(thisStruct);
        logTable = [logTable;thisRow];
    end

end


%% Stack to string
function stackText = xASL_test_StackToString(stack)

    % Fallback
    stackText = '';
    
    % Iterate over elements
    for iElement = 1:size(stack,1)
        thisFile = stack(iElement).file;
        thisName = stack(iElement).name;
        thisLine = stack(iElement).line;
        stackText = [stackText ', ' thisName ': line ' num2str(thisLine)];
    end
    
    % Remove initial ' ,'
    if ~isempty(stackText) && length(stackText)>3
        stackText = stackText(3:end);
    end

end


%% Write logging information to XLSX table
function xASL_test_writeLoggingTable(filePath,logTable)

    % Write file
    try
        writetable(logTable,filePath,'WriteVariableNames',true);
    catch ME
        fprintf(2,'Writing log.xlsx failed: %s\n', ME.message);
    end

end



