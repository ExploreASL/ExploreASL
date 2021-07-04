function UnitTest = xASL_ut_UnitTest_function_PSNR(TestRepository)
%xASL_ut_UnitTest_function_PSNR Individual unit test for xASL_stat_PSNR
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
% EXAMPLE:      UnitTests(1) = xASL_ut_UnitTest_function_PSNR(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL

%% Initialize test structure

% Insert test name here
UnitTest.name = 'xASL_stat_PSNR';

% Define whether you are testing a module, submodule or function
UnitTest.unit = 'Function';

%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Artificial noise examples';

% Start the test
testTime = tic;

% Prepare test
xASL_Copy(fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_2_0','asl\001_asl.nii.gz'),...
          fullfile(TestRepository,'UnitTesting','working_directory','test.nii.gz'))
im = xASL_io_Nifti2Im(fullfile(TestRepository,'UnitTesting','working_directory','test.nii.gz'));
xASL_delete(fullfile(TestRepository,'UnitTesting','working_directory','test.nii.gz'));

% Run your test here
imNoiseless = im;
randomMatrix = rand(size(im),'single');
downscaledMatrix = randomMatrix.*0.00001;
noisyImage1 = imNoiseless+downscaledMatrix*1;
noisyImage2 = imNoiseless+downscaledMatrix*10;
noisyImage3 = imNoiseless+downscaledMatrix*100;
noisyImage4 = imNoiseless+downscaledMatrix*1000;
PSNR1 = xASL_stat_PSNR(imNoiseless,noisyImage1);
PSNR2 = xASL_stat_PSNR(imNoiseless,noisyImage2);
PSNR3 = xASL_stat_PSNR(imNoiseless,noisyImage3);
PSNR4 = xASL_stat_PSNR(imNoiseless,noisyImage4);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~isnumeric(PSNR1) || ~isnumeric(PSNR2) || ~isnumeric(PSNR3) || ~isnumeric(PSNR4)
    testCondition = false;
end
if ~(PSNR1>PSNR2) || ~(PSNR2>PSNR3) || ~(PSNR3>PSNR4)
    testCondition = false;
end

clear PSNR1 PSNR2 PSNR3 PSNR4 im imNoiseless noisyImage1 noisyImage2 noisyImage3 noisyImage4 randomMatrix downscaledMatrix

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% Test run 2

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Dummy examples';

% Start the test
testTime = tic;

% Prepare test
xASL_Copy(fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_2_0','asl\001_asl.nii.gz'),...
          fullfile(TestRepository,'UnitTesting','working_directory','test.nii.gz'))
im = xASL_io_Nifti2Im(fullfile(TestRepository,'UnitTesting','working_directory','test.nii.gz'));
xASL_delete(fullfile(TestRepository,'UnitTesting','working_directory','test.nii.gz'));

% Run your test here
PSNR1 = xASL_stat_PSNR(0,0);
PSNR2 = xASL_stat_PSNR(im,im);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~isnan(PSNR1)
    testCondition = false;
end
if ~isinf(PSNR2)
    testCondition = false;
end

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


