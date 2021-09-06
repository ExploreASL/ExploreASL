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
y1 = xASL_stat_StdNan(A);
y2 = xASL_stat_StdNan(B);
y3 = xASL_stat_StdNan(C);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~(mean(y1)-stdA<0.001)
    testCondition = false;
end
if ~(mean(mean(y2))-stdB<0.001)
    testCondition = false;
end
if ~(mean(mean(mean(y3)))-stdC<0.001)
    testCondition = false;
end

% Clean-Up
clear A B C stdA stdB stdC y1 y2 y3

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;


%% Test run 3

% Give your individual subtest a name
UnitTest.tests(3).testname = 'Multidimensional test with NaN values (2D)';

% Start the test
testTime = tic;

% Run your test here
A = zeros(300,100);
B = ones(300,100);
C = NaN(300,100);
Z = [A,B,C];
Z = Z'; % => Z(:,1) = (0,...,0,1,...,1,NaN,...,NaN)
y = xASL_stat_StdNan(Z);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~((mean(y)-0.5)<0.01)
    testCondition = false;
end

% Get test duration
UnitTest.tests(3).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(3).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


