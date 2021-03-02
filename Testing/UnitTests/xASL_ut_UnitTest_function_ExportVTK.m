function UnitTest = xASL_ut_UnitTest_function_ExportVTK(TestRepository)
%xASL_ut_UnitTest_function_ExportVTK Individual unit test for xASL_io_ExportVTK
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
% EXAMPLE:      UnitTests(1) = xASL_ut_UnitTest_function_ExportVTK(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL

%% Initialize test structure

% Insert test name here
UnitTest.name = 'xASL_io_ExportVTK';

% Define whether you are testing a module, submodule or function
UnitTest.unit = 'Function';

%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'NIFTI path given';

% Start the test
testTime = tic;

% Run your test here
workingDirectory = fullfile(TestRepository,'UnitTesting','working_directory');

% Initialization
[x] = ExploreASL_Initialize([],0);

% Test NIFTIs
testNiftiVTK = fullfile(x.MyPath,'External','TestDataSet','analysis','Sub-001','ASL_1','M0.nii.gz');

% Run test
xASL_Copy(testNiftiVTK,fullfile(workingDirectory,'image.nii'),true);
xASL_io_ExportVTK(x.MyPath,fullfile(workingDirectory,'image.nii'));

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~exist(fullfile(workingDirectory,'image.nii'),'file')
    testCondition = false;
end
if ~exist(fullfile(workingDirectory,'export.vtk'),'file')
    testCondition = false;
end

% Delete test files
xASL_delete(fullfile(workingDirectory,'image.nii'));
xASL_delete(fullfile(workingDirectory,'export.vtk'));

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% Test run 2

% Give your individual subtest a name
UnitTest.tests(2).testname = 'NIFTI image & export path given';

% Start the test
testTime = tic;

% Run your test here
workingDirectory = fullfile(TestRepository,'UnitTesting','working_directory');

% Initialization
[x] = ExploreASL_Initialize([],0);

% Test NIFTIs
testNiftiVTK = fullfile(x.MyPath,'External','TestDataSet','analysis','Sub-001','ASL_1','M0.nii.gz');

% Run test
xASL_Copy(testNiftiVTK,fullfile(workingDirectory,'image.nii'),true);
image = xASL_io_Nifti2Im(fullfile(workingDirectory,'image.nii'));
xASL_io_ExportVTK(x.MyPath,image,[],fullfile(workingDirectory,'image.vtk'));

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~exist(fullfile(workingDirectory,'image.nii'),'file')
    testCondition = false;
end
if ~exist(fullfile(workingDirectory,'image.vtk'),'file')
    testCondition = false;
end

% Delete test files
xASL_delete(fullfile(workingDirectory,'image.nii'));
xASL_delete(fullfile(workingDirectory,'image.vtk'));

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;


%% Test run 3

% Give your individual subtest a name
UnitTest.tests(3).testname = 'Neither NIFTI path nor export path given';

% Start the test
testTime = tic;

% Run your test here
workingDirectory = fullfile(TestRepository,'UnitTesting','working_directory');

% Initialization
[x] = ExploreASL_Initialize([],0);

% Test NIFTIs
testNiftiVTK = fullfile(x.MyPath,'External','TestDataSet','analysis','Sub-001','ASL_1','M0.nii.gz');

% Run test
xASL_Copy(testNiftiVTK,fullfile(workingDirectory,'image.nii'),true);
image = xASL_io_Nifti2Im(fullfile(workingDirectory,'image.nii'));
xASL_io_ExportVTK(x.MyPath,image);

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~exist(fullfile(workingDirectory,'image.nii'),'file')
    testCondition = false;
end
if ~exist(fullfile(pwd,'export.vtk'),'file')
    testCondition = false;
end

% Delete test files
xASL_delete(fullfile(workingDirectory,'image.nii'));
xASL_delete(fullfile(pwd,'export.vtk'));

% Get test duration
UnitTest.tests(3).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(3).passed = testCondition;

%% Test run 4

% Give your individual subtest a name
UnitTest.tests(4).testname = 'Masking & Resampling';

% Start the test
testTime = tic;

% Run your test here
workingDirectory = fullfile(TestRepository,'UnitTesting','working_directory');

% Initialization
[x] = ExploreASL_Initialize([],0);

% Test NIFTIs
testNiftiVTK = fullfile(x.MyPath,'External','TestDataSet','analysis','Sub-001','ASL_1','M0.nii.gz');
testMask = fullfile(x.MyPath,'External','SPMmodified','toolbox','cat12','templates_volumes','brainmask.nii');

% Run test
xASL_Copy(testNiftiVTK,fullfile(workingDirectory,'image.nii'),true);
xASL_io_ExportVTK(x.MyPath,fullfile(workingDirectory,'image.nii'),testMask);

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~exist(fullfile(workingDirectory,'image.nii'),'file')
    testCondition = false;
end
if ~exist(fullfile(workingDirectory,'export.vtk'),'file')
    testCondition = false;
end

% Delete test files
xASL_delete(fullfile(workingDirectory,'image.nii'));
xASL_delete(fullfile(workingDirectory,'export.vtk'));

% Get test duration
UnitTest.tests(4).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(4).passed = testCondition;


%% Test run 5

% Give your individual subtest a name
UnitTest.tests(5).testname = 'No input arguments';

% Start the test
testTime = tic;

% Run test
[x] = ExploreASL_Initialize([],0);
xASL_io_ExportVTK(x.MyPath);
[warnMessage, ~] = lastwarn();

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~strcmp(warnMessage,'Neither input image path nor matrix given...')
    testCondition = false;
end

% Get test duration
UnitTest.tests(5).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(5).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end





