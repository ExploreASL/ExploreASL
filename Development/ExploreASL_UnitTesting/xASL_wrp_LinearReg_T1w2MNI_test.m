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

%% Test 1: xASL_wrp_LinearReg_T1w2MNI: TestDataSet (low quality)

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
ProcessData = true;
iWorker = 1;
nWorkers = 1;

% Initialize the data set
x = ExploreASL_Initialize(DataParPath, ProcessData, iWorker, nWorkers);

% Add necessary additional x structure fields
x.SUBJECTDIR      = '<ROOT>/<SUBJECT>';
x.LockDir         = ['<ROOT>/lock/xASL_module_Structural/<SUBJECT>'];

curFolder = fullfile(testDir,'TestFolder','TestDataSet','Sub-001');
x.P.Path_FLAIR = fullfile(curFolder,'FLAIR.nii');
x.P.Path_WMH_SEGM = fullfile(curFolder,'WMH_SEGM.nii');
x.P.Path_T1 = fullfile(curFolder,'T1.nii');
x.P.Path_ASL4D = fullfile(curFolder,'ASL_1\ASL4D.nii');
x.P.Path_M0 = fullfile(curFolder,'ASL_1\M0.nii');
x.P.Path_ASL4D_RevPE = fullfile(curFolder,'ASL_1\ASL4D_RevPE.nii');

% Unzip Nifti files
IM = xASL_io_Nifti2Im(x.P.Path_T1);
if size(IM,4)>1 || size(IM,5)>1 || size(IM,6)>1 || size(IM,7)>1
    warning(['Too many dims, using first: ' x.P.Path_T1]);
    xASL_Copy(x.P.Path_T1, x.P.Path_T1_ORI, true);
    xASL_io_SaveNifti(x.P.Path_T1, x.P.Path_T1, IM(:,:,:,1,1,1,1), [], false);
end

IM = xASL_io_Nifti2Im(x.P.Path_FLAIR);
if size(IM,4)>1 || size(IM,5)>1 || size(IM,6)>1 || size(IM,7)>1
    warning(['Too many dims, using first: ' x.P.Path_T1]);
    xASL_Move(x.P.Path_FLAIR, x.P.Path_FLAIR_ORI, true);
    xASL_io_SaveNifti(x.P.Path_FLAIR, x.P.Path_FLAIR, IM(:,:,:,1,1,1,1), [], false);
end

% Run test
xASL_wrp_LinearReg_T1w2MNI(x);

% Check quality setting
assert(exist(fullfile(testDir,'TestFolder','TestDataSet','Sub-001','ASL_1','ASL4D.mat'),'file')==2)

% What could be tested here? (WORK IN PROGRESS)





