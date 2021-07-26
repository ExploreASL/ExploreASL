function UnitTest = xASL_ut_function_xASL_adm_DeleteFileList(TestRepository)
%xASL_ut_function_xASL_adm_DeleteFileList Individual unit test for xASL_adm_DeleteFileList
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_adm_DeleteFileList(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Delete file lists (paths with spaces)';

% Start the test
testTime = tic;

% Prepare tests

% Create six test files in test directory
xASL_adm_CreateDir(fullfile(TestRepository,'UnitTesting','working_directory','dir with spaces'));
fid01 = fopen(fullfile(TestRepository,'UnitTesting','working_directory','dir with spaces','test_1.txt'), 'wt');
fid02 = fopen(fullfile(TestRepository,'UnitTesting','working_directory','dir with spaces','test_2.txt'), 'wt');
fid03 = fopen(fullfile(TestRepository,'UnitTesting','working_directory','dir with spaces','test_3.txt'), 'wt');
fid04 = fopen(fullfile(TestRepository,'UnitTesting','working_directory','dir with spaces','dont_1.json'), 'wt');
fid05 = fopen(fullfile(TestRepository,'UnitTesting','working_directory','dir with spaces','dont_2.json'), 'wt');
fid06 = fopen(fullfile(TestRepository,'UnitTesting','working_directory','dir with spaces','dont_3.json'), 'wt');
fclose('all');

% Run your test here
pathWithSpaces = fullfile(TestRepository,'UnitTesting','working_directory','dir with spaces');
[deletedFilesA,NotDeletedFilesA] = xASL_adm_DeleteFileList(pathWithSpaces, '^.+.txt$');

rmdir(fullfile(TestRepository,'UnitTesting','working_directory','dir with spaces'), 's');

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if isempty(deletedFilesA)
    testCondition = false;
end
if ~isempty(NotDeletedFilesA)
    testCondition = false; 
end
if ~(numel(deletedFilesA)==3)
    testCondition = false;
end
if ~(numel(NotDeletedFilesA)==0)
    testCondition = false;
end

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% Test run 2

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Delete file lists (paths without spaces)';

% Start the test
testTime = tic;

% Prepare tests

% Create six test files in test directory
xASL_adm_CreateDir(fullfile(TestRepository,'UnitTesting','working_directory','dir_without_spaces'));
fid07 = fopen(fullfile(TestRepository,'UnitTesting','working_directory','dir_without_spaces','test_1.txt'), 'wt');
fid08 = fopen(fullfile(TestRepository,'UnitTesting','working_directory','dir_without_spaces','test_2.txt'), 'wt');
fid09 = fopen(fullfile(TestRepository,'UnitTesting','working_directory','dir_without_spaces','test_3.txt'), 'wt');
fid10 = fopen(fullfile(TestRepository,'UnitTesting','working_directory','dir_without_spaces','dont_1.json'), 'wt');
fid11 = fopen(fullfile(TestRepository,'UnitTesting','working_directory','dir_without_spaces','dont_2.json'), 'wt');
fid12 = fopen(fullfile(TestRepository,'UnitTesting','working_directory','dir_without_spaces','dont_3.json'), 'wt');
fclose('all');


% Run your test here
pathWithoutSpaces = fullfile(TestRepository,'UnitTesting','working_directory','dir_without_spaces');
[deletedFilesB,NotDeletedFilesB] = xASL_adm_DeleteFileList(pathWithoutSpaces, '^.+.txt$');

rmdir(fullfile(TestRepository,'UnitTesting','working_directory','dir_without_spaces'), 's');

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if isempty(deletedFilesB)
    testCondition = false;
end
if ~isempty(NotDeletedFilesB)
    testCondition = false; 
end
if ~(numel(deletedFilesB)==3)
    testCondition = false;
end
if ~(numel(NotDeletedFilesB)==0)
    testCondition = false;
end

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;




%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end

