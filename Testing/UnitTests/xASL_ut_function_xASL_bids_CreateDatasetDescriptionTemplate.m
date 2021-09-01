function UnitTest = xASL_ut_function_xASL_bids_CreateDatasetDescriptionTemplate(TestRepository)
%xASL_ut_function_xASL_bids_CreateDatasetDescriptionTemplate Individual unit test for xASL_bids_CreateDatasetDescriptionTemplate
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_bids_CreateDatasetDescriptionTemplate(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Empty draft';

% Start the test
testTime = tic;

% Run your test here
draft = struct;
testVersion = '1.2.3';
[json] = xASL_bids_CreateDatasetDescriptionTemplate(draft,testVersion);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~isstruct(draft)
    testCondition = false;
end
% Check required fields
if ~isfield(json,'Name') || ...
   ~isfield(json,'BIDSVersion')
    testCondition = false;
end

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% Test run 2

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Draft with fields that are not allowed';

% Start the test
testTime = tic;

% Run your test here
draft = struct;
draft.NotAllowed = 'Test';
testVersion = '1.2.3';
[json] = xASL_bids_CreateDatasetDescriptionTemplate(draft,testVersion);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~isstruct(draft)
    testCondition = false;
end
% Check required fields
if ~isfield(json,'Name') || ...
   ~isfield(json,'BIDSVersion')
    testCondition = false;
end
% Make sure that not allowed field was removed
if isfield(json,'NotAllowed')
    testCondition = false;
end

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


