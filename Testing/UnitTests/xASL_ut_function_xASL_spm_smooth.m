function UnitTest = xASL_ut_function_xASL_spm_smooth(TestRepository)
%xASL_ut_function_xASL_spm_smooth Individual unit test for xASL_spm_smooth
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_spm_smooth(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Smooth a T1w image using different kernel sizes';

% Start the test
testTime = tic;

% Set-up test NIfTI
testDirsAndFiles.testSource = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_3_0_pre_release','sub-001','anat');
testDirsAndFiles.testDestination = fullfile(TestRepository,'UnitTesting','working_directory','testSubject');
testDirsAndFiles.testFile = fullfile(testDirsAndFiles.testDestination,'test.nii.gz');
testDirsAndFiles.testFileResult1 = fullfile(testDirsAndFiles.testDestination,'result1.nii.gz');
testDirsAndFiles.testFileResult2 = fullfile(testDirsAndFiles.testDestination,'result2.nii.gz');
testDirsAndFiles.testFileResult3 = fullfile(testDirsAndFiles.testDestination,'result3.nii.gz');

% Set up directory and copy NIfTI
xASL_adm_CreateDir(testDirsAndFiles.testDestination);
xASL_Copy(fullfile(testDirsAndFiles.testSource,'sub-001_acq-003_T1w.nii.gz'),testDirsAndFiles.testFile,1);

% Fallback
testCondition = true;

% Test: Smooth
try
    xASL_spm_smooth(testDirsAndFiles.testFile,[4,4,4],testDirsAndFiles.testFileResult1);
    xASL_spm_smooth(testDirsAndFiles.testFile,[8,8,8],testDirsAndFiles.testFileResult2);
    xASL_spm_smooth(testDirsAndFiles.testFile,[16,16,16],testDirsAndFiles.testFileResult3);
catch ME
    warning(ME.identifier, '%s', ME.message);
    testCondition = false;
    diary off;
end

% Test: Load data
try
    imOriginal = xASL_io_Nifti2Im(testDirsAndFiles.testFile);
    im1 = xASL_io_Nifti2Im(testDirsAndFiles.testFileResult1);
    im2 = xASL_io_Nifti2Im(testDirsAndFiles.testFileResult2);
    im3 = xASL_io_Nifti2Im(testDirsAndFiles.testFileResult3);
catch ME
    imOriginal = [];
    im1 = [];
    im2 = [];
    im3 = [];
    warning(ME.identifier, '%s', ME.message);
    testCondition = false;
    diary off;
end

% Define one or multiple test conditions here
if isempty(imOriginal) || isempty(im1) || isempty(im2) || isempty(im3)
    testCondition = false;
end

% Check image size first
if isequal(size(imOriginal),size(im1)) && isequal(size(imOriginal),size(im2)) && isequal(size(imOriginal),size(im3))

    % The sum of image values compared to the original image should stay almost the same
    if abs(sum(im1(:))/sum(imOriginal(:))-1) > 0.01
        testCondition = false;
    end
    if abs(sum(im2(:))/sum(imOriginal(:))-1) > 0.01
        testCondition = false;
    end
    if abs(sum(im3(:))/sum(imOriginal(:))-1) > 0.01
        testCondition = false;
    end

    % The RMSE should be pretty small as well
    if xASL_ut_GetRMSE(imOriginal, im1) > 1.0
        testCondition = false;
    end
    if xASL_ut_GetRMSE(imOriginal, im2) > 1.0
        testCondition = false;
    end
    if xASL_ut_GetRMSE(imOriginal, im3) > 1.0
        testCondition = false;
    end
    
    % I also expect that the variance decreases for gaussian blurring / smoothing
    resultValues.varOriginal = var(imOriginal(:));
    resultValues.varIm1 = var(im1(:));
    resultValues.varIm2 = var(im2(:));
    resultValues.varIm3 = var(im3(:));
    
    if ~(resultValues.varOriginal>resultValues.varIm1) ...
        || ~(resultValues.varIm1>resultValues.varIm2) ...
        || ~(resultValues.varIm2>resultValues.varIm3) 
        testCondition = false;
    end
else
    testCondition = false;
end



% Delete test data
xASL_delete(testDirsAndFiles.testDestination,true);

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

