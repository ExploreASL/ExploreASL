function UnitTest = xASL_ut_function_xASL_adm_CheckFileCount(TestRepository)
%xASL_ut_function_xASL_adm_CheckFileCount Individual unit test for xASL_adm_CheckFileCount
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_adm_CheckFileCount(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Find three files';

% Start the test
testTime = tic;

% Create three test files in test directory
fid1 = fopen(fullfile(TestRepository,'UnitTesting','working_directory','test_1.txt'), 'wt');
fid2 = fopen(fullfile(TestRepository,'UnitTesting','working_directory','test_2.txt'), 'wt');
fid3 = fopen(fullfile(TestRepository,'UnitTesting','working_directory','test_3.txt'), 'wt');

% Run your test here
[result, files] = xASL_adm_CheckFileCount(fullfile(TestRepository,'UnitTesting','working_directory'),'test*',3,0);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~result
    testCondition = false;
end
if ~(numel(files)==3)
    testCondition = false; 
end

% Delete test files
xASL_delete(fullfile(TestRepository,'UnitTesting','working_directory','test_1.txt'))
xASL_delete(fullfile(TestRepository,'UnitTesting','working_directory','test_2.txt'))
xASL_delete(fullfile(TestRepository,'UnitTesting','working_directory','test_3.txt'))

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% Test run 2

% Give your individual subtest a name
UnitTest.tests(2).testname = 'No matching files';

% Start the test
testTime = tic;

% Run your test here
[result, files] = xASL_adm_CheckFileCount(fullfile(TestRepository,'UnitTesting','working_directory'),'test*',3,0);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if result
    testCondition = false;
end
if ~(numel(files)==0)
    testCondition = false; 
end

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;



%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


