function UnitTest = xASL_ut_UnitTest_module_Structural(TestRepository)
%xASL_ut_UnitTest_module_Structural Individual unit test for xASL_module_Structural
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
% EXAMPLE:      UnitTests(1) = xASL_ut_UnitTest_module_Structural(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL

%% Initialize test structure

% Insert test name here
UnitTest.name = 'xASL_module_Structural';

% Define whether you are testing a module, submodule or function
UnitTest.unit = 'Module';

%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Default settings';

% Start the test
testTime = tic;

% Prepare the test data
pathTestPatient = fullfile(TestRepository,'UnitTesting','working_directory','Patient');
xASL_Copy(fullfile(TestRepository,'UnitTesting','dro_files','0_01_Derivatives_Data','Patient'), pathTestPatient);

% Run your test here
x = ExploreASL(pathTestPatient,0,[1 0 0]);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~isstruct(x)
    testCondition = false;
end

% Add more conditions here ... (compare data with reference patient dataset 1_00_Structural_Module)

% Delete test data
xASL_delete(pathTestPatient,true)

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


