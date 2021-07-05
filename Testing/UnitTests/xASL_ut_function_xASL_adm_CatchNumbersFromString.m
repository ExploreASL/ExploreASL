function UnitTest = xASL_ut_function_xASL_adm_CatchNumbersFromString(TestRepository)
%xASL_ut_UnitTest_function_CatchNumbersFromString Individual unit test for xASL_adm_CatchNumbersFromString
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_adm_CatchNumbersFromString(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Extract integer example';

% Start the test
testTime = tic;

% Run your test here
[number] = xASL_adm_CatchNumbersFromString('abc123def');

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~isnumeric(number)
    testCondition = false;
end
if isnumeric(number) && ~number==123
    testCondition = false; 
end

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% Test run 2

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Extract double example';

% Start the test
testTime = tic;

% Run your test here
[number] = xASL_adm_CatchNumbersFromString('abc123.456def');

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~isnumeric(number)
    testCondition = false;
end
if isnumeric(number) && ~(number(1)==123 && number(2)==456)
    testCondition = false; 
end

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;


%% Test run 3

% Give your individual subtest a name
UnitTest.tests(3).testname = 'No number example';

% Start the test
testTime = tic;

% Run your test here
[number] = xASL_adm_CatchNumbersFromString('abcdef');

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~isnumeric(number)
    testCondition = false;
end
if isnumeric(number) && ~isnan(number)
    testCondition = false;
end

% Get test duration
UnitTest.tests(3).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(3).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


