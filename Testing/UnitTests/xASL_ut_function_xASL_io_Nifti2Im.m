function UnitTest = xASL_ut_function_xASL_io_Nifti2Im(TestRepository)
%xASL_ut_function_xASL_io_Nifti2Im Individual unit test for xASL_Nifti2Im
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_io_Nifti2Im(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Load DRO image matrix (default options)';

% Start the test
testTime = tic;

% Run your test here
aslFile = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_3_0','rawdata','sub-001','perf','sub-001_acq-001_asl.nii.gz');
testFile = fullfile(TestRepository,'UnitTesting','working_directory','001_asl.nii.gz');
unzippedTestFile = fullfile(TestRepository,'UnitTesting','working_directory','001_asl.nii');

% Copy test NIFTI
copyfile(aslFile,testFile);

% Run function
imTest = xASL_io_Nifti2Im(testFile);

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~exist(unzippedTestFile,'file')
    testCondition = false; % Test failed
end
if ~isa(imTest,'single')
    testCondition = false; % Test failed
end
if numel(imTest)~=(64*64*12*3) % Matrix size: [64, 64, 12, 3]
    testCondition = false; % Test failed
end

% Remove file after test
delete(unzippedTestFile);

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


