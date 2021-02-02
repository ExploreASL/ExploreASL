function UnitTest = xASL_ut_UnitTest_Template(TestRepository)
%xASL_ut_UnitTest_Template Individual unit test template
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
% DESCRIPTION:  To improve quality control, each developer can define unit
%               tests using this template. Please save the individual unit
%               test in the "Tests" directory.
%               Only modify the UnitTest.name, UnitTest.unit and the
%               individual UnitTest.tests(x).testname fields.
%               Insert your test code/functions in the "Run your test here ..."
%               areas.
%               Naming conventions: please name your individual unit test
%               according to the following scheme:
%               xASL_qc_UnitTest_[name of the module/submodule/function]
%               For example: the unit test of the xASL_module_ASL would be called:
%               xASL_ut_UnitTest_module_ASL
%
% EXAMPLE:      UnitTests(1) = xASL_ut_UnitTest_Template(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL

%% Initialize test structure

% Insert test name here
UnitTest.name = 'Template Test';

% Define whether you are testing a module, submodule or function
UnitTest.unit = 'Module';

%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Check A';

% Start the test
testTime = tic;

% Run your test here
% ...

% Define one or multiple test conditions here
testCondition = true;

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;




%% Test run 2

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Check B';

% Start the test
testTime = tic;

% Run your test here
% ...

% Define one or multiple test conditions here
testCondition = true;

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;

%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end





