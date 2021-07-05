function UnitTest = xASL_ut_function_xASL_adm_CompareLists(TestRepository)
%xASL_ut_function_xASL_adm_CompareLists Individual unit test for xASL_adm_CompareLists
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_adm_CompareLists(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Simple list example';

% Start the test
testTime = tic;

% Run your test here
testListA = {1, 2, 3, 'A', 'B', 'C'}';
testListB = {4, 5, 6, 'A', 'B', 'C'}';
[newList] = xASL_adm_CompareLists(testListA, testListB);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~isempty(newList{1,1}) || ~isempty(newList{2,1}) || ~isempty(newList{3,1}) || ...
   ~isempty(newList{1,2}) || ~isempty(newList{2,2}) || ~isempty(newList{3,2})
    testCondition = false;
end
if ~strcmp(newList{4,1},'A') || ~strcmp(newList{5,1},'B') || ~strcmp(newList{6,1},'C') || ...
   ~newList{4,2} || ~newList{5,2} || ~newList{6,2}
end

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% Test run 2

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Different list length example';

% Start the test
testTime = tic;

% Run your test here
testListA = {1, 2, 3, 'A', 'B', 'C'}';
testListB = {4, 5, 6, 'A', 'B'}';
[newList] = xASL_adm_CompareLists(testListA, testListB);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~isempty(newList{1,1}) || ~isempty(newList{2,1}) || ~isempty(newList{3,1}) || ...
   ~isempty(newList{1,2}) || ~isempty(newList{2,2}) || ~isempty(newList{3,2})
    testCondition = false;
end
if ~strcmp(newList{4,1},'A') || ~strcmp(newList{5,1},'B') || ~strcmp(newList{6,1},'C') || ...
   ~newList{4,2} || ~newList{5,2} || newList{6,2}
end

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


