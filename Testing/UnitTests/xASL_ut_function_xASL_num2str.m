function UnitTest = xASL_ut_function_xASL_num2str(TestRepository)
%xASL_ut_function_xASL_num2str Individual unit test for xASL_num2str
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_str2num(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Row vector, column vector, scalar number';

% Start the test
testTime = tic;

% Run your test here
outA = xASL_num2str([1, 2, 3]');
outB = xASL_num2str([1, 2, 3]);
outC = xASL_num2str(1);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~ischar(outA) || ~ischar(outB) || ~ischar(outC)
    testCondition = false;
end
if isempty(regexp(outA,'1,2,3', 'once'))
    testCondition = false;
end
if isempty(regexp(outB,'1,2,3', 'once'))
    testCondition = false;
end
if ~strcmp(outC,'1')
    testCondition = false;
end

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;



%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


