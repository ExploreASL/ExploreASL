function UnitTest = xASL_ut_function_xASL_im_CompareNiftis(TestRepository)
%xASL_ut_function_xASL_im_CompareNiftis Individual unit test for xASL_im_CompareNiftis
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_im_CompareNiftis(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Identical images';

% Start the test
testTime = tic;

% Run your test here
aslFile = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_3_0','rawdata','sub-001','perf','sub-001_acq-001_asl.nii.gz');
testFile = fullfile(TestRepository,'UnitTesting','working_directory','001_asl.nii.gz');
testFile2 = fullfile(TestRepository,'UnitTesting','working_directory','002_asl.nii.gz');

% Copy files
xASL_Copy(aslFile,testFile);
xASL_Copy(aslFile,testFile2);

% Run function
[identical,RMSE,minDiff,maxDiff,dimCheck] = xASL_im_CompareNiftis(testFile,testFile2,1);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~identical
    testCondition = false;
end
if RMSE>0.001
    testCondition = false;
end
if minDiff>0.001
    testCondition = false;
end
if maxDiff>0.001
    testCondition = false;
end
if ~dimCheck
    testCondition = false;
end

% Remove file after test
xASL_delete(testFile);
xASL_delete(testFile2);

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% Test run 2

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Different images';

% Start the test
testTime = tic;

% Run your test here
aslFile = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_3_0','rawdata','sub-001','ground_truth','sub-001_acq-002_Perfmap.nii.gz');
aslFile2 = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_3_0','rawdata','sub-001','ground_truth','sub-001_acq-002_T1map.nii.gz');
testFile = fullfile(TestRepository,'UnitTesting','working_directory','perf.nii.gz');
testFile2 = fullfile(TestRepository,'UnitTesting','working_directory','t1.nii.gz');

% Copy files
xASL_Copy(aslFile,testFile);
xASL_Copy(aslFile2,testFile2);

% Run function
[identical,RMSE,minDiff,maxDiff,dimCheck] = xASL_im_CompareNiftis(testFile,testFile2,1);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if identical
    testCondition = false;
end
if RMSE<0.001
    testCondition = false;
end
if ~dimCheck
    testCondition = false;
end

% Remove file after test
xASL_delete(testFile);
xASL_delete(testFile2);

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


