%xASL_wrp_LinearReg_T1w2MNI_test Script to test the xASL_wrp_LinearReg_T1w2MNI function
%
% FORMAT:       RESULT = runtests('xASL_wrp_LinearReg_T1w2MNI_test');
% 
% INPUT:        None
%
% OUTPUT:       Console window
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This script is supposed to run various unit tests of xASL_wrp_LinearReg_T1w2MNI:
%
%           1) Run a test using the default TestDataSet inputs with low quality setting
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLES: RESULT = runtests('xASL_wrp_LinearReg_T1w2MNI_test');
% __________________________________
% Copyright 2015-2020 ExploreASL

% Check if parameter file exists
test_parameter_file = fullfile(pwd,'Development','ExploreASL_UnitTesting','xASL_test_parameters.json');
if xASL_exist(test_parameter_file)
    val = jsondecode(fileread(test_parameter_file));
    if isfield(val,'xASLdir'), xASLdir = val.xASLdir; end
    if isfield(val,'testDir'), testDir = val.testDir; end
    if strcmp(xASLdir,'pwd')
        xASLdir = pwd;
    end
end

% PRECONDITIONS
TestDataSetPath = fullfile(testDir,'TestFolder','TestDataSet');

%% Test 1: TestDataSet (low quality)

% Check if folder was already created
if exist(fullfile(testDir,'TestFolder'),'dir')==7
    fprintf('Remove existing folder %s...\n', testDir)
    rmdir(fullfile(testDir,'TestFolder'),'s');
end
fprintf('Creating test folder in %s...\n', testDir)
mkdir(fullfile(testDir,'TestFolder'))

xASL_Copy(fullfile(xASLdir,'External\TestDataSet'), fullfile(testDir,'TestFolder','TestDataSet'))
fprintf('Copy test data to %s...\n', fullfile(testDir,'TestFolder','TestDataSet'))

% Initialize test
fprintf('Initialize test input...\n')
DataParPath = fullfile(TestDataSetPath,'DataParameters_LowQ.json');

% Initialize the data set
x = ExploreASL_Initialize(DataParPath, true, 1, 1);
x = xASL_qc_Default_Test_Initialize(x,TestDataSetPath,x.SUBJECTS{1});
x = xASL_init_FileSystem(x);

% Unzip T1 and FLAIR
xASL_qc_Default_Test_Unzip(x);

% Backup
xBackup = x;

% Run test
xASL_wrp_LinearReg_T1w2MNI(x);

% What could be tested here? (WORK IN PROGRESS)

[match, er1, er2] = comp_struct(x,xBackup,1,0,1e-3);

% assert()




