function UnitTest = xASL_ut_function_xASL_test_GetLogContent(TestRepository)
%xASL_ut_function_xASL_test_GetLogContent Individual unit test for xASL_test_GetLogContent
%
% INPUT:        TestRepository - Path to test repository.
%
% OUTPUT:       UnitTest  - Test structure
%               name      - Name of tested module or submodule (char array)
%               unit      - Insert one of the following: 'Module', 'Submodule' or 'Function'
%               passed    - Result of all subtests combined (true or false)
%               test      - Structure with individual subtest results
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Should be run using xASL_ut_UnitTesting.
%
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_test_GetLogContent(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Read in log files (no ouput file)';

% Start the test
testTime = tic;

% Run your test here
testDirectory = fullfile(TestRepository,'UnitTesting','io_files');

% Read test files
[logContent] = xASL_test_GetLogContent(testDirectory,0,1,0);

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~istable(logContent)
    testCondition = false; % Test failed
end
if ~ischar(logContent.Module{1})
    testCondition = false; % Test failed
end
if ~ischar(logContent.Subject{1})
    testCondition = false; % Test failed
end
if ~ischar(logContent.Message{1})
    testCondition = false; % Test failed
end
if ~ischar(logContent.File{1})
    testCondition = false; % Test failed
end
if ~ischar(logContent.Line{1})
    testCondition = false; % Test failed
end

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;



%% Test run 2

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Read in log files (ouput tsv file)';

% Start the test
testTime = tic;

% Run your test here
testDirectory = fullfile(TestRepository,'UnitTesting','io_files');

% Output file
outputFile = fullfile(TestRepository,'UnitTesting','io_files','logContent.tsv');

% Read test files
[logContent] = xASL_test_GetLogContent(testDirectory,0,1,1);

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~istable(logContent)
    testCondition = false; % Test failed
end
if ~ischar(logContent.Module{1})
    testCondition = false; % Test failed
end
if ~ischar(logContent.Subject{1})
    testCondition = false; % Test failed
end
if ~ischar(logContent.Message{1})
    testCondition = false; % Test failed
end
if ~ischar(logContent.File{1})
    testCondition = false; % Test failed
end
if ~ischar(logContent.Line{1})
    testCondition = false; % Test failed
end
if ~exist(outputFile,'file')
    testCondition = false; % Test failed
end

% Remove file after test
delete(outputFile);

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;



%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


