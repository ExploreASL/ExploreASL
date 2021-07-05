function UnitTest = xASL_ut_function_xASL_tsvRead(TestRepository)
%xASL_ut_function_xASL_tsvRead Individual unit test for xASL_tsvRead
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_tsvRead(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Read test file (default options)';

% Start the test
testTime = tic;

% Run your test here
testFile = fullfile(TestRepository,'UnitTesting','io_files','TestFile.tsv');
CellContents = xASL_tsvRead(testFile);

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~iscell(CellContents)
	testCondition = false; % Test failed
end
if ~(size(CellContents,1)==5 && size(CellContents,2)==3)
    testCondition = false; % Test failed
end
if ~(strcmp(CellContents{1,1},'Parameter') && strcmp(CellContents{1,2},'Description') && strcmp(CellContents{1,3},'Value'))
    testCondition = false; % Test failed
end
if ~(ischar(CellContents{2,1}) && ischar(CellContents{2,2}) && isa(CellContents{2,3},'double'))
    testCondition = false; % Test failed
end

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% Test run 2

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Read test file (bStruct option)';

% Start the test
testTime = tic;

% Run your test here
testFile = fullfile(TestRepository,'UnitTesting','io_files','TestFile.tsv');
CellContents = xASL_tsvRead(testFile,true);

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~isstruct(CellContents)
	testCondition = false; % Test failed
end
if ~(numel(CellContents.Parameter)==4 && numel(CellContents.Description)==4)
    testCondition = false; % Test failed
end
if ~(strcmp(CellContents.Parameter{1,1},'Height') && strcmp(CellContents.Parameter{2,1},'Width') && strcmp(CellContents.Parameter{3,1},'Length') && strcmp(CellContents.Parameter{4,1},'Depth'))
    testCondition = false; % Test failed
end
if ~(ischar(CellContents.Parameter{1,1}) && ischar(CellContents.Description{1,1}) && isa(CellContents.Value(1),'double'))
    testCondition = false; % Test failed
end

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


