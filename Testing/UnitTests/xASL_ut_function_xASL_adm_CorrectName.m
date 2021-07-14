function UnitTest = xASL_ut_function_xASL_adm_CorrectName(TestRepository)
%xASL_ut_function_xASL_adm_CorrectName Individual unit test for xASL_adm_CorrectName
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_adm_CorrectName(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Test these symbols: $%&';

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


%% Test run 2

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Test these symbols: []()?!$%&';

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

%% Test run 3

% Give your individual subtest a name
UnitTest.tests(3).testname = 'Test these symbols: ÄÖÜäöü/\';

% Start the test
testTime = tic;

% Run your test here
strOut = xASL_adm_CorrectName('abcÄÖÜäöü/\/\/\def');

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~strcmp(strOut,'abc_def')
    testCondition = false; 
end

% Get test duration
UnitTest.tests(3).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(3).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


