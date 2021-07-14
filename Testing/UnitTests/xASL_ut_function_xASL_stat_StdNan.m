function UnitTest = xASL_ut_function_xASL_stat_StdNan(TestRepository)
%xASL_ut_function_xASL_stat_StdNan Individual unit test for xASL_adm_CheckFileCount
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_stat_StdNan(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Multiple test calls with different arrays (1D)';

% Start the test
testTime = tic;

% Run your test here
y1 = xASL_stat_StdNan([0, 0, 0, 0, 0]);
y2 = xASL_stat_StdNan([0, 1, 2, 3, 4]);
y3 = xASL_stat_StdNan([0, NaN, 0, NaN, 0]);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~(y1==0)
    testCondition = false;
end
if ~(y2==std([0, 1, 2, 3, 4]))
    testCondition = false;
end
if ~(y3==0)
    testCondition = false;
end

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% Test run 2

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Multiple test calls with different arrays (2D, 3D, 4D)';

% Start the test
testTime = tic;

% Run your test here
A = randi([0 1],100,100,'double');
B = randi([0 1],100,100,100,'double');
C = randi([0 1],100,100,100,100,'double');
stdA = std(A(:));
stdB = std(B(:));
stdC = std(C(:));
nanA = NaN(100,100);
nanB = NaN(100,100,100);
nanC = NaN(100,100,100,100);
concatA = [A, nanA];
concatB = [B, nanB];
concatC = [C, nanC];
y1 = xASL_stat_StdNan(concatA);
y2 = xASL_stat_StdNan(concatB);
y3 = xASL_stat_StdNan(concatC);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~(y1==stdA)
    testCondition = false;
end
if ~(y2==stdB)
    testCondition = false;
end
if ~(y3==stdC)
    testCondition = false;
end

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


