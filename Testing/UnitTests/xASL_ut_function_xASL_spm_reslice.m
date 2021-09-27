function UnitTest = xASL_ut_function_xASL_spm_reslice(TestRepository)
%xASL_ut_function_xASL_spm_reslice Individual unit test for xASL_spm_reslice
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_spm_reslice(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Reslice a T1w image ...';

% Start the test
testTime = tic;

% Set-up test NIfTI
testDirsAndFiles.testSource = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_3_0_pre_release','sub-001','anat');
testDirsAndFiles.testDestination = fullfile(TestRepository,'UnitTesting','working_directory','testSubject');
testDirsAndFiles.testFile = fullfile(testDirsAndFiles.testDestination,'test.nii.gz');
testDirsAndFiles.testFileResult1 = fullfile(testDirsAndFiles.testDestination,'result1.nii.gz');

% Set up directory and copy NIfTI
xASL_adm_CreateDir(testDirsAndFiles.testDestination);
xASL_Copy(fullfile(testDirsAndFiles.testSource,'sub-001_acq-003_T1w.nii.gz'),testDirsAndFiles.testFile,1);

% Fallback
testCondition = true;

% Test: Reslice
try
    % xASL_spm_reslice(testDirsAndFiles.testFile,...);
catch ME
    warning(ME.identifier, '%s', ME.message);
    testCondition = false;
    diary off;
end

% Test: Load data
try
    im1 = xASL_io_Nifti2Im(testDirsAndFiles.testFileResult1);
catch ME
    im1 = [];
    warning(ME.identifier, '%s', ME.message);
    testCondition = false;
    diary off;
end

% Define one or multiple test conditions here




% Delete test data
xASL_delete(testDirsAndFiles.testDestination),true)

% Clean-up
clearvars -except UnitTest TestRepository testCondition testTime 

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


%% Determine RMSE of two images
function RMSE = xASL_ut_GetRMSE(imageA, imageB)
    
    if isequal(size(imageA),size(imageB))
        % Calculate
        RMSE = sqrt(mean((imageA(:) - imageB(:)).^2))*2/sqrt(mean(abs(imageA(:)) + abs(imageB(:))).^2);
    else
        % Resample A to B
        [imageA] = xASL_im_ResampleIM(imageA, [1 0 0 0;0 1 0 0;0 0 1 0; 0 0 0 1], [1 0 0 0;0 1 0 0;0 0 1 0; 0 0 0 1], size(imageB));
        % Calculate
        RMSE = sqrt(mean((imageA(:) - imageB(:)).^2))*2/sqrt(mean(abs(imageA(:)) + abs(imageB(:))).^2);
    end
end

