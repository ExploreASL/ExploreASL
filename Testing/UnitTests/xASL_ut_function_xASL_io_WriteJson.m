function UnitTest = xASL_ut_function_xASL_io_WriteJson(TestRepository)
%xASL_ut_function_xASL_io_WriteJson Individual unit test for xASL_io_WriteJson
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_io_WriteJson(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2023 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Write a test file';

% Start the test
testTime = tic;

% Run your test here
testFile = fullfile(TestRepository, 'UnitTesting', 'io_files', 'testFile.json');
testFileWrite = fullfile(TestRepository, 'UnitTesting', 'working_directory', 'testFile.json');

% Read and write the test file
jsonRef = xASL_io_ReadJson(testFile);
xASL_io_WriteJson(testFileWrite, jsonRef);

% Read the written file again and compare with the original file
jsonTest = xASL_io_ReadJson(testFileWrite);

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~exist(testFileWrite,'file')
    testCondition = false; % Test failed
end
if ~isstruct(jsonTest)
	testCondition = false; % Test failed
end
if ~isequal(jsonRef, jsonTest)    
	testCondition = false; % Test failed
end

% Remove file after test
delete(testFileWrite);

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% Test run 2

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Read and write test file (force overwrite)';

% Start the test
testTime = tic;

% Run your test here
testFile = fullfile(TestRepository, 'UnitTesting', 'io_files', 'testFile.json');
testFileWrite = fullfile(TestRepository, 'UnitTesting', 'working_directory', 'testFile.json');

% Copy a JSON file to the output
xASL_Copy(testFile, testFileWrite);

jsonRef = struct();
jsonRef.a = 1;
jsonRef.b = double([2.1, 3.1])';
jsonRef.c = 'test';

% Write test file while overwriting the current content
xASL_io_WriteJson(testFileWrite, jsonRef);

% Read test file
jsonTest = xASL_io_ReadJson(testFileWrite);

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~exist(testFileWrite,'file')
    testCondition = false; % Test failed
end
if ~isstruct(jsonTest)
	testCondition = false; % Test failed
end
if ~isequal(jsonRef, jsonTest)    
	testCondition = false; % Test failed
end

% Remove file after test
delete(testFileWrite);

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end
