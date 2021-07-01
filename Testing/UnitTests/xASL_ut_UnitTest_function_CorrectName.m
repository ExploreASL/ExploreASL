function UnitTest = xASL_ut_UnitTest_function_CorrectName(TestRepository)
%xASL_ut_UnitTest_function_CorrectName Individual unit test for xASL_adm_CorrectName
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
% EXAMPLE:      UnitTests(1) = xASL_ut_UnitTest_function_CorrectName(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL

%% Initialize test structure

% Insert test name here
UnitTest.name = 'xASL_adm_CatchNumbersFromString';

% Define whether you are testing a module, submodule or function
UnitTest.unit = 'Function';

%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Simple text example';

% Start the test
testTime = tic;

% Run your test here
strOut = xASL_adm_CorrectName('abc$%&def');

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~isvarname(strOut)
    testCondition = false;
end
if isvarname(strOut) && ~strcmp(strOut,'abc_def')
    testCondition = false; 
end

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Another text example';

% Start the test
testTime = tic;

% Run your test here
strOut = xASL_adm_CorrectName('[]()???!!!abc$%&def');

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~strcmp(strOut,'_abc_def')
    testCondition = false; 
end

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


