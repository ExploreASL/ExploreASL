function [UnitTests,UnitTestsTable] = xASL_test_UnitTesting(bPull)
%xASL_test_UnitTesting Main script to run all individual unit tests
%
% FORMAT: [UnitTests,UnitTestsTable] = xASL_test_UnitTesting([bPull])
%
% INPUT:        bPull           - Pull up-to-date testing repository (BOOLEAN, OPTIONAL, DEFAULT = true)
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
% EXAMPLE:      [UnitTests,UnitTestsTable] = xASL_test_UnitTesting;
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL

    %% Admin
	if nargin < 1 || isempty(bPull)
		bPull = true;
	end
	
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
    clearvars -except xASL_dir bPull

    %% Get Test Repository
    
    % Try to find the ExploreASL/Testing repository on the same folder level
    [scriptPath,~,~] = fileparts(mfilename('fullpath'));
    potentialTestingDirectory = strrep(scriptPath,fullfile('ExploreASL','Testing'),'Testing');
    if exist(potentialTestingDirectory,'dir')
        fprintf('Testing repository found...\n');
        TestRepository = potentialTestingDirectory;
    else
        TestRepository = [];
    end
    
    % Give user feedback if repository was not found
    if isempty(TestRepository)
        fprintf(2,'ExploreASL was unable to find the Testing repository in your local directory...\n');
        UnitTests = NaN;
        UnitTestsTable = NaN;
        return
    end
	
	% Update test repository
	cd(TestRepository); % Go to test repository
	if bPull
		system('git fetch','-echo');
		system('git pull','-echo');
	end
	cd(xASL_dir) % Go back to ExploreASL
    
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
    
    % Get working directory for unit tests
    workingDirectory = fullfile(TestRepository,'UnitTesting','working_directory');
    
    % Unit test template
    UnitTestTemplate.name = 'Template';
    UnitTestTemplate.unit = 'Function';
    UnitTestTemplate.tests = struct;
    UnitTestTemplate.passed = true;
    
    % Get unit tests
    addpath(fullfile(scriptPath,'UnitTests'));
    fileList = dir(fullfile(scriptPath,'UnitTests','*.m'));
    
    % Iterate over tests
    for test = 1:size(fileList,1)
        % Make sure the working directory is empty
        xASL_delete(workingDirectory,true);
        xASL_adm_CreateDir(workingDirectory);
        % Get test handle
        testScript = fileList(test).name;
        [~,testScript,~] = fileparts(testScript);
        testHandle = str2func(testScript);
        % Run test
        try
            UnitTest = testHandle(TestRepository);
        catch
            UnitTest = struct;
            UnitTest.tests = struct;
            UnitTest.passed = false;
        end
        % Assign unit test name and unit definitions
        UnitTest = xASL_ut_InitUnitTest(UnitTest,testHandle);
        % Write test to array
        UnitTests(test) = UnitTest;
    end
    
    % Order struct fields
    UnitTests = orderfields(UnitTests, UnitTestTemplate);
    
    %% Export TSV
    for iTest = 1:size(UnitTests,2)
        UnitTestsCells{iTest,1} = UnitTests(iTest).name;
        UnitTestsCells{iTest,2} = UnitTests(iTest).unit;
        UnitTestsCells{iTest,3} = num2str(UnitTests(iTest).passed);
    end
    xASL_tsvWrite(UnitTestsCells,fullfile(TestRepository,'results.tsv'),1);
    
    %% Export table as well
    UnitTestsTable = struct2table(UnitTests);
    
    %% Print test results
    clc
    fprintf('====================================== TEST RESULTS ======================================\n\n')
    disp(UnitTestsTable);
    fprintf('==========================================================================================\n')

end


%% Initialize unit test
function UnitTest = xASL_ut_InitUnitTest(UnitTest,testHandle)

    % Get name of script
    testObject = functions(testHandle);

    % Insert test name here
    UnitTest.name = testObject.function(regexp(testObject.function,'_xASL_')+1:end);

    % Define whether you are testing a module, submodule or function
    UnitTest.unit = testObject.function(regexp(testObject.function,'xASL_ut_')+8:regexp(testObject.function,'_xASL_')-1);

end



