function UnitTest = xASL_ut_UnitTest_function_QuantileNan(TestRepository)
%xASL_ut_UnitTest_function_QuantileNan Individual unit test for xASL_adm_CheckFileCount
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
% EXAMPLE:      UnitTests(1) = xASL_ut_UnitTest_function_QuantileNan(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL

%% Initialize test structure

% Insert test name here
UnitTest.name = 'xASL_stat_QuantileNan';

% Define whether you are testing a module, submodule or function
UnitTest.unit = 'Function';

%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Three basic examples';

% Start the test
testTime = tic;

% Run your test here
x1 = [0, 0, 0, 0, 0, 0, 0];
x2 = [1, 2, 3, 4, 5, 6, 7];
x3 = [NaN, 1, NaN, 2, NaN, 3, NaN, 4, NaN, 5, NaN, 6, NaN, 7];
y1 = xASL_stat_QuantileNan(x1,0.5);
y2 = xASL_stat_QuantileNan(x2,0.5);
y3 = xASL_stat_QuantileNan(x3,0.5);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~(y1==0)
    testCondition = false;
end
if ~(y2==4) % Reference from wolfram alpha
    testCondition = false;
end
if ~(y3==4) % Reference from wolfram alpha
    testCondition = false;
end

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


