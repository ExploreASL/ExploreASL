function UnitTest = xASL_ut_function_xASL_imp_StudyParPriority(TestRepository)
%xASL_ut_function_xASL_imp Individual unit test for xASL_imp_StudyParPriority
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_imp(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________



%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Read test file (default options)';

% Start the test
testTime = tic;

% Run your test here
testFile = fullfile(TestRepository,'UnitTesting','io_files','studyParMultiParameter.json');
testJson = xASL_io_ReadDataPar(testFile);

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~isfield(testJson,'StudyPars')
	testCondition = false; % Test failed
end
if length(testJson.StudyPars) ~= 4
    testCondition = false; % Test failed
end

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;

%% Test run 2 - run 3

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Test studyParameters of run 3';

% Start the test
testTime = tic;

% Run your test here
studyPar = xASL_imp_StudyParPriority(testJson, 'gamma', '1', '3');

% Define one or multiple test conditions here
testCondition = true; % Fallback
if studyPar.PostLabelingDelay ~= 2
	testCondition = false; % Test failed
end

if studyPar.BackgroundSuppression ~= false
	testCondition = false; % Test failed
end

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;

%% Test run 3 a non-special studyParameters

% Give your individual subtest a name
UnitTest.tests(3).testname = 'Test a general non-special study parameters';

% Start the test
testTime = tic;

% Run your test here
studyPar = xASL_imp_StudyParPriority(testJson, 'gamma', '1', '1');

% Define one or multiple test conditions here
testCondition = true; % Fallback
if studyPar.PostLabelingDelay ~= 2
	testCondition = false; % Test failed
end

if studyPar.BackgroundSuppression ~= true
	testCondition = false; % Test failed
end

% Get test duration
UnitTest.tests(3).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(3).passed = testCondition;

%% Test run 4 - a non-defined study parameters

% Give your individual subtest a name
UnitTest.tests(4).testname = 'Test non-restricted study parameters';

% Start the test
testTime = tic;

% Run your test here
studyPar = xASL_imp_StudyParPriority(testJson, '', '', '');

% Define one or multiple test conditions here
testCondition = true; % Fallback
if studyPar.PostLabelingDelay ~= 4
	testCondition = false; % Test failed
end

if studyPar.BackgroundSuppression ~= false
	testCondition = false; % Test failed
end

% Get test duration
UnitTest.tests(4).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(4).passed = testCondition;

%% Test run 5

% Give your individual subtest a name
UnitTest.tests(5).testname = 'Test general study parameters';

% Start the test
testTime = tic;

% Run your test here
studyPar = xASL_imp_StudyParPriority(testJson, 'alpha2', '1', '3');

% Define one or multiple test conditions here
testCondition = true; % Fallback
if studyPar.PostLabelingDelay ~= 3
	testCondition = false; % Test failed
end

if studyPar.BackgroundSuppression ~= false
	testCondition = false; % Test failed
end

% Get test duration
UnitTest.tests(5).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(5).passed = testCondition;

%% Test run 6

% Give your individual subtest a name
UnitTest.tests(6).testname = 'Test session/run study parameters';

% Start the test
testTime = tic;

% Run your test here
studyPar = xASL_imp_StudyParPriority(testJson, '', '2', '1');

% Define one or multiple test conditions here
testCondition = true; % Fallback
if studyPar.PostLabelingDelay ~= 4
	testCondition = false; % Test failed
end

if studyPar.BackgroundSuppression ~= true
	testCondition = false; % Test failed
end

% Get test duration
UnitTest.tests(6).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(6).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


