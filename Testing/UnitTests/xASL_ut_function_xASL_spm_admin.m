function UnitTest = xASL_ut_function_xASL_spm_admin(TestRepository)
%xASL_ut_function_xASL_spm_admin Individual unit test for xASL_spm_admin
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_spm_admin(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Read gzipped NIfTI using SPM admin';

% Start the test
testTime = tic;

% Set-up test NIfTI
testDirsAndFiles.testSource = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_3_0','rawdata','sub-001','anat');
testDirsAndFiles.testDestination = fullfile(TestRepository,'UnitTesting','working_directory','testSubject');
testDirsAndFiles.testFile = fullfile(testDirsAndFiles.testDestination,'test.nii.gz');

% Set up directory and copy NIfTI
xASL_adm_CreateDir(testDirsAndFiles.testDestination);
xASL_Copy(fullfile(testDirsAndFiles.testSource,'sub-001_acq-003_T1w.nii.gz'),testDirsAndFiles.testFile,1);

% Fallback
testCondition = true;

% Test: Admin
try
    pathOut = xASL_spm_admin(testDirsAndFiles.testFile);
catch ME
    warning(ME.identifier, '%s', ME.message);
    testCondition = false;
    diary off;
end

% Test: Load data
try
    imTest = xASL_io_Nifti2Im(testDirsAndFiles.testFile);
catch ME
    imTest = [];
    warning(ME.identifier, '%s', ME.message);
    testCondition = false;
    diary off;
end
fileList = xASL_adm_GetFileList(testDirsAndFiles.testDestination);

% Define one or multiple test conditions here
if isempty(imTest)
    testCondition = false;
end
if isempty(fileList)
    testCondition = false;
else
    [result.dir, result.name, result.extension] = fileparts(fileList{1});
    % NIfTI should be unzipped now
    if ~strcmp(result.extension,'.nii')
        testCondition = false;
    end
    % NIfTI should be called test now
    if ~strcmp(result.name,'test')
        testCondition = false;
    end
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



