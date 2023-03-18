function UnitTest = xASL_ut_function_xASL_io_ReadJson(TestRepository)
%xASL_ut_function_xASL_io_ReadJson Individual unit test for xASL_io_ReadJson
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_io_ReadJson(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2023 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Read test file';

% Start the test
testTime = tic;

% Run your test here
testFile = fullfile(TestRepository, 'UnitTesting', 'io_files', 'testFile.json');
CellContents = xASL_io_ReadJson(testFile);

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~(isstruct(CellContents))
    testCondition = false; % Test failed
end
if ~(isfield(CellContents, 'folderHierarchy') && isfield(CellContents, 'tokenOrdering') && isfield(CellContents, 'tokenSessionAliases') && isfield(CellContents, 'tokenScanAliases') && isfield(CellContents, 'bMatchDirectories'))
    testCondition = false; % Test failed
end
if ~(iscell(CellContents.folderHierarchy) && isnumeric(CellContents.tokenOrdering) && iscell(CellContents.tokenSessionAliases))
    testCondition = false; % Test failed
end
if ~(size(CellContents.tokenSessionAliases,1) == 4 && size(CellContents.tokenOrdering,1) == 4)
    testCondition = false; % Test failed
end

if ~(isequal(strcmp(CellContents.folderHierarchy{2}, '^(ASL|T1w|M0|T2|FLAIR)(\d)$'), 1))
	testCondition = false; % Test failed
end

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;

%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end
