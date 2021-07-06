function UnitTest = xASL_ut_function_xASL_tsvWrite(TestRepository)
%xASL_ut_function_xASL_tsvWrite Individual unit test for xASL_tsvWrite
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_tsvWrite(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Read and write test file (default options)';

% Start the test
testTime = tic;

% Run your test here
testFile = fullfile(TestRepository,'UnitTesting','io_files','TestFile.tsv');
testFileWrite = fullfile(TestRepository,'UnitTesting','working_directory','TestFile.tsv');
InputCell = xASL_tsvRead(testFile);
% Write test file
xASL_tsvWrite(InputCell, testFileWrite, false, false);
% Read test file
OutputCell = xASL_tsvRead(testFileWrite);

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~exist(testFileWrite,'file')
    testCondition = false; % Test failed
end
if ~iscell(OutputCell)
	testCondition = false; % Test failed
end
if ~(size(OutputCell,1)==5 && size(OutputCell,2)==3)
    testCondition = false; % Test failed
end
if ~(strcmp(OutputCell{1,1},'Parameter') && strcmp(OutputCell{1,2},'Description') && strcmp(OutputCell{1,3},'Value'))
    testCondition = false; % Test failed
end
if ~(ischar(OutputCell{2,1}) && ischar(OutputCell{2,2}) && isa(OutputCell{2,3},'double'))
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
testFile = fullfile(TestRepository,'UnitTesting','io_files','TestFile.tsv');
testFileWrite = fullfile(TestRepository,'UnitTesting','working_directory','TestFile.tsv');
InputCell = xASL_tsvRead(testFile);
% Write test file
xASL_tsvWrite(InputCell, testFileWrite, false, false);
xASL_tsvWrite(InputCell, testFileWrite, true, false); % Overwrite file from previous step
% Read test file
OutputCell = xASL_tsvRead(testFileWrite);

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~exist(testFileWrite,'file')
    testCondition = false; % Test failed
end
if ~iscell(OutputCell)
	testCondition = false; % Test failed
end
if ~(size(OutputCell,1)==5 && size(OutputCell,2)==3)
    testCondition = false; % Test failed
end
if ~(strcmp(OutputCell{1,1},'Parameter') && strcmp(OutputCell{1,2},'Description') && strcmp(OutputCell{1,3},'Value'))
    testCondition = false; % Test failed
end
if ~(ischar(OutputCell{2,1}) && ischar(OutputCell{2,2}) && isa(OutputCell{2,3},'double'))
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


