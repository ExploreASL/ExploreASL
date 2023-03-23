function UnitTest = xASL_ut_function_xASL_adm_SortStringNumbers(TestRepository)
%xASL_ut_function_xASL_adm_SortStringNumbers Individual unit test for xASL_io_WriteJson
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_adm_SortStringNumbers(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2023 ExploreASL

%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Test string number sorting';

% Start the test
testTime = tic;

% Run your test here
NewList{1,1} = 'ASL_1';
NewList{2,1} = 'ASL_6';
NewList{3,1} = 'ASL_3';
NewList{4,1} = 'ASL_05';
NewList{5,1} = 'ASL_12';

% Sort the list
SortedList = xASL_adm_SortStringNumbers(NewList);

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~strcmp(SortedList{1},'ASL_1')
    testCondition = false; % Test failed
end
if ~strcmp(SortedList{2},'ASL_3')
    testCondition = false; % Test failed
end
if ~strcmp(SortedList{3},'ASL_05')
    testCondition = false; % Test failed
end
if ~strcmp(SortedList{4},'ASL_6')
    testCondition = false; % Test failed
end
if ~strcmp(SortedList{5},'ASL_12')
    testCondition = false; % Test failed
end

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;

%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end
