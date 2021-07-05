function UnitTest = xASL_ut_function_xASL_str2num(TestRepository)
%xASL_ut_function_xASL_str2num Individual unit test for xASL_str2num
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
UnitTest.tests(1).testname = 'Number types: integer, double, positive and negative';

% Start the test
testTime = tic;

% Run your test here
numOutA = xASL_str2num('123');
numOutB = xASL_str2num('-123');
numOutC = xASL_str2num('123.456');
numOutD = xASL_str2num('-123.456');

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~isnumeric(numOutA) || ~isnumeric(numOutB) || ~isnumeric(numOutC) || ~isnumeric(numOutD)
    testCondition = false;
end
if isnumeric(numOutA) && isnumeric(numOutB) && isnumeric(numOutC) && isnumeric(numOutD)
    if ~(numOutA==123) || ~(numOutB==-123) || ~(numOutC==123.456) || ~(numOutD==-123.456)
        testCondition = false;
    end
end

% Clean up
clear numOutA numOutB numOutC numOutD

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;

%% Test run 1

% Give your individual subtest a name
UnitTest.tests(2).testname = 'List of numbers';

% Start the test
testTime = tic;

% Run your test here
testArr = xASL_str2num({'1','2','3'},0);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~isnumeric(testArr)
    testCondition = false;
end
if isnumeric(testArr) && numel(testArr)==3
    if ~(testArr(1)==1) || ~(testArr(2)==2) || ~(testArr(3)==3)
        testCondition = false;
    end
end

% Clean up
clear testArr

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


