function [UnitTests,UnitTestsTable] = xASL_ut_UnitTesting
%xASL_qc_UnitTesting Main script to run all individual unit tests
%
% INPUT:        n/a
%
% OUTPUT:       UnitTests       - structure containing the unit test results
%               UnitTestsTable  - table containing the unit test results
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function is used to run all individual unit tests. To
%               define a unit test, please use the xASL_qc_UnitTest_Template.
%               The idea is that this script can run independently from the
%               rest of ExploreASL, to enable unbiased and robust testing.
% 				This means you, as a developer, should not use/add ExploreASL
% 				functions within this script!
%
% EXAMPLE:      [UnitTests,UnitTestsTable] = xASL_ut_UnitTesting;
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

    %% Clean Up
    clc
    clear

    %% Get Test Repository
    
    % Try to find the ExploreASL/Testing repository on the same folder level
    [scriptPath,~,~] = fileparts(mfilename('fullpath'));
    potentialTestingDirectory = strrep(scriptPath,fullfile('ExploreASL','Testing'),'Testing');
    if isfolder(potentialTestingDirectory)
        fprintf('Testing repository found...\n');
        TestRepository = potentialTestingDirectory;
    else    
        TestRepository = [];
    end
    
    if isempty(TestRepository)
        if usejava('desktop')
            TestRepository = uigetdir([],'Select test repository...');
        else
            TestRepository = input('Insert test repository: ');
        end
        if ~ischar(TestRepository)
            error('Test repository path required...');
        end
    end

    %% Test Workflow

    % Get unit tests
    addpath(fullfile(scriptPath,'UnitTests'));
    fileList = dir(fullfile(scriptPath,'UnitTests','*.m'));
    
    % Iterate over tests
    for test = 1:size(fileList,1)
        testScript = fileList(test).name;
        [~,testScript,~] = fileparts(testScript);
        testHandle = str2func(testScript);
        UnitTests(test) = testHandle(TestRepository);
    end
        
    
    
    
    %% Export table as well
    UnitTestsTable = struct2table(UnitTests);
    
    %% Print test results
    clc
    fprintf('================================= TEST RESULTS =================================\n\n')
    disp(UnitTestsTable);
    fprintf('================================================================================\n')

end



