function UnitTest = xASL_ut_function_xASL_im_ResampleLinear(TestRepository)
%xASL_ut_function_xASL_im_ResampleLinear Individual unit test for xASL_im_ResampleLinear
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_im_ResampleLinear(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL

%% Paths

paths.sourcePath = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_3_0','rawdata','sub-001','anat','sub-001_acq-003_T1w.nii.gz');
paths.testPath = fullfile(TestRepository,'UnitTesting','working_directory','test.nii.gz');
paths.testPathUnzipped = fullfile(TestRepository,'UnitTesting','working_directory','test.nii');

%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Resample 3D';

% Start the test
testTime = tic;

% Run your test here
xASL_Copy(paths.sourcePath,paths.testPath);
image = xASL_io_Nifti2Im(paths.testPath);
output = xASL_im_ResampleLinear(image,[100,200,300]);
xASL_delete(paths.testPathUnzipped);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~isnumeric(output)
    testCondition = false;
end
if ~isequal(size(output),[100,200,300])
    testCondition = false; 
end

% Clean up
clear image output

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% Test run 2

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Resample 2D';

% Start the test
testTime = tic;

% Run your test here
xASL_Copy(paths.sourcePath,paths.testPath);
image = xASL_io_Nifti2Im(paths.testPath);
image = image(:,:,1);
output = xASL_im_ResampleLinear(image,[100,200]);
xASL_delete(paths.testPathUnzipped);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~isnumeric(output)
    testCondition = false;
end
if ~isequal(size(output),[100,200])
    testCondition = false; 
end

% Clean up
clear image output

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;

%% Test run 3

% Give your individual subtest a name
UnitTest.tests(3).testname = 'Resample 1D';

% Start the test
testTime = tic;

% Run your test here
xASL_Copy(paths.sourcePath,paths.testPath);
image = xASL_io_Nifti2Im(paths.testPath);
image = image(:,1,1);
output = xASL_im_ResampleLinear(image,[100]);
xASL_delete(paths.testPathUnzipped);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~isnumeric(output)
    testCondition = false;
end
if numel(output)~=100
    testCondition = false;
end

% Clean up
clear image output

% Get test duration
UnitTest.tests(3).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(3).passed = testCondition;

%% Test run 4

% Give your individual subtest a name
UnitTest.tests(4).testname = 'Dimension mismatch';

% Start the test
testTime = tic;

% Run your test here
xASL_Copy(paths.sourcePath,paths.testPath);
image = xASL_io_Nifti2Im(paths.testPath);
image = image(:,:,1);
message = '';
try
    output = xASL_im_ResampleLinear(image,[100,200,300]);
catch ME
    message = ME.message;
end
xASL_delete(paths.testPathUnzipped);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if isempty(message)
    testCondition = false;
end
if ~strcmp(message,'Mismatch of new size and image dimension...')
    testCondition = false;
end

% Clean up
clear image output

% Get test duration
UnitTest.tests(4).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(4).passed = testCondition;


%% Test run 5

% Give your individual subtest a name
UnitTest.tests(5).testname = 'Resample 3D (up, down, etc.)';

% Start the test
testTime = tic;

% Run your test here
imageA = zeros(100,25,100);
imageB = ones(100,25,100).*1;
imageC = ones(100,25,100).*2;
imageD = ones(100,25,100).*3;
image = [imageA imageB imageC imageD];
outputUp = xASL_im_ResampleLinear(image,[200,200,200]);
outputDown = xASL_im_ResampleLinear(image,[50,50,50]);

% Edge cases
output1x1x1 = xASL_im_ResampleLinear(image,[1,1,1]);
output2x2x2 = xASL_im_ResampleLinear(image,[2,2,2]);
output100x100x1 = xASL_im_ResampleLinear(image,[100,100,1]);
output1x100x1 = xASL_im_ResampleLinear(image,[1,100,1]);

% Mean
meanSource = mean(image,'all');
meanUp = mean(outputUp,'all');
meanDown = mean(outputDown,'all');
stdSource = std(image(:)');
stdUp = std(outputUp(:)');
stdDown = std(outputDown(:)');

% Tolerance
tolMean = 0.00001;
tolStd = 0.02;

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~isnumeric(outputUp) || ~isnumeric(outputDown)
    testCondition = false;
end
if ~isequal(size(outputUp),[200,200,200]) || ~isequal(size(outputDown),[50,50,50])
    testCondition = false; 
end
if abs(meanSource-meanUp)>tolMean || abs(meanSource-meanDown)>tolMean
    testCondition = false; 
end
if abs(stdSource-stdUp)>tolStd || abs(stdSource-stdDown)>tolStd
    testCondition = false; 
end

% Edge cases
if ~isnumeric(output1x1x1) || isnan(output1x1x1)
    testCondition = false;
end
if ~isnumeric(output2x2x2)
    testCondition = false;
end
if ~isnumeric(output100x100x1) || sum(isnan(output100x100x1(:)'))>0
    testCondition = false;
end
if ~isnumeric(output1x100x1) || sum(isnan(output1x100x1(:)'))>0
    testCondition = false;
end

% Clean up
clear image outputUp outputDown

% Get test duration
UnitTest.tests(5).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(5).passed = testCondition;



%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


