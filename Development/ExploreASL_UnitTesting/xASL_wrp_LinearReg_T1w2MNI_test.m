%xASL_wrp_LinearReg_T1w2MNI_test Script to test the xASL_wrp_LinearReg_T1w2MNI function
%
% FORMAT:       RESULT = runtests('xASL_wrp_LinearReg_T1w2MNI_test');
% 
% INPUT:        None
%
% OUTPUT:       Console window
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This script is supposed to run various unit tests of xASL_Initialize:
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
DataParPath = fullfile(testDir,'TestFolder','TestDataSet','DataParameters_LowQ.json');

% Initialize the data set
x = ExploreASL_Initialize(DataParPath, true, 1, 1);
x = xASL_qc_Default_Test_Initialize(x);

% Add paths which should actually be defined in: x = xASL_init_FileSystem(x);
curFolder = fullfile(testDir,'TestFolder','TestDataSet','Sub-001');
x.P.Path_FLAIR = fullfile(curFolder,'FLAIR.nii');
x.P.Path_WMH_SEGM = fullfile(curFolder,'WMH_SEGM.nii');
x.P.Path_T1 = fullfile(curFolder,'T1.nii');
x.P.Path_ASL4D = fullfile(curFolder,'ASL_1\ASL4D.nii');
x.P.Path_M0 = fullfile(curFolder,'ASL_1\M0.nii');
x.P.Path_ASL4D_RevPE = fullfile(curFolder,'ASL_1\ASL4D_RevPE.nii');

% Unzip T1 and FLAIR
xASL_qc_Default_Test_Unzip(x);

% Run test
xASL_wrp_LinearReg_T1w2MNI(x);

% Check quality setting
assert(exist(fullfile(testDir,'TestFolder','TestDataSet','Sub-001','ASL_1','ASL4D.mat'),'file')==2)

% What could be tested here? (WORK IN PROGRESS)





