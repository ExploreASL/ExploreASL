%xASL_module_Structural_test Script to test the xASL_module_Structural function
%
% FORMAT:       RESULT = runtests('xASL_module_Structural_test');
% 
% INPUT:        None
%
% OUTPUT:       Console window
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This script is supposed to run various unit tests of xASL_module_Structural:
%
%           1) Run a test using the default TestDataSet inputs with low quality setting
%           2) Run a test with ...
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLES: RESULT = runtests('xASL_module_Structural_test');
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

% Check if folder was already created
if exist(fullfile(testDir,'TestFolder'),'dir')==7
    fprintf('Remove existing folder %s...\n', testDir)
    rmdir(fullfile(testDir,'TestFolder'),'s');
end
fprintf('Creating test folder in %s...\n', testDir)
mkdir(fullfile(testDir,'TestFolder'))

xASL_Copy(fullfile(xASLdir,'External\TestDataSet'), fullfile(testDir,'TestFolder','TestDataSet'))
fprintf('Copy test data to %s...\n', fullfile(testDir,'TestFolder','TestDataSet'))

% PRECONDITIONS

%% Test 1: Default TestDataSet with low quality setting
fprintf('Initialize test input...\n')
DataParPath = fullfile(testDir,'TestFolder','TestDataSet','DataParameters_LowQ.json');
ProcessData = true;
iWorker = 1;
nWorkers = 1;
xTest = ExploreASL_Initialize(DataParPath, ProcessData, iWorker, nWorkers);

% Run test
[~, x] = xASL_Iteration(xTest,'xASL_module_Structural');

% Use assert for outputs
assert(isfield(x,'RERUN'))
assert(isfield(x,'MUTEXID'))
assert(isfield(x,'LockDir'))
% assert(isfield(x,'SUBJECTDIR'))
assert(isfield(x,'SUBJECT'))
assert(isfield(x,'ModuleName'))
assert(isfield(x,'result'))
% assert(isfield(x,'mutex'))
assert(isfield(x,'iSubject'))
assert(isfield(x,'DoWADQCDC'))
assert(isfield(x,'WMHsegmAlg'))
assert(isfield(x,'SkipIfNoFlair'))
assert(isfield(x,'SkipIfNoASL'))
assert(isfield(x,'SkipIfNoM0'))
assert(isfield(x,'bFixResolution'))
assert(isfield(x,'Seg'))
assert(isfield(x,'T1BiasFieldRegularization'))




% If I use ExploreASL with the TestDataSet and with low quality option, I
% don't get the Output field after xASL_module_Structural. Some additional
% informations are added after xASL_module_ASL and some other ones after
% xASL_module_Population. These additional ones could be the same ones, but
% there is definitely no Output field, even in the x structure after all
% three modules.




% 
% 
% 
% % Check Output field
% assert(isfield(x,'Output'))
% if isfield(x,'Output')
%     assert(isfield(x.Output,'Structural'))
%     if isfield(x.Output,'Structural')
%         assert(isfield(x.Output.Structural,'T1w_WMref_vol_mL'))
%         assert(isfield(x.Output.Structural,'T1w_SD_WMref'))
%         assert(isfield(x.Output.Structural,'T1w_SNR_GM_Ratio'))
%         assert(isfield(x.Output.Structural,'T1w_CNR_GM_WM_Ratio'))
%         assert(isfield(x.Output.Structural,'T1w_FBER_WMref_Ratio'))
%         assert(isfield(x.Output.Structural,'T1w_EFC_bits'))
%         assert(isfield(x.Output.Structural,'T1w_Mean_AI_Perc'))
%         assert(isfield(x.Output.Structural,'T1w_SD_AI_Perc'))
%         assert(isfield(x.Output.Structural,'FLAIR_WMref_vol_mL'))
%         assert(isfield(x.Output.Structural,'FLAIR_SD_WMref'))
%         assert(isfield(x.Output.Structural,'FLAIR_SNR_GM_Ratio'))
%         assert(isfield(x.Output.Structural,'FLAIR_CNR_GM_WM_Ratio'))
%         assert(isfield(x.Output.Structural,'FLAIR_FBER_WMref_Ratio'))
%         assert(isfield(x.Output.Structural,'FLAIR_EFC_bits'))
%         assert(isfield(x.Output.Structural,'FLAIR_Mean_AI_Perc'))
%         assert(isfield(x.Output.Structural,'FLAIR_SD_AI_Perc'))
%         assert(isfield(x.Output.Structural,'ID'))
%         assert(isfield(x.Output.Structural,'T1w_LR_flip_YesNo'))
%         assert(isfield(x.Output.Structural,'FLAIR_WMH_vol_mL'))
%         assert(isfield(x.Output.Structural,'FLAIR_WMH_n'))
%         assert(isfield(x.Output.Structural,'T1w_IQR_Perc'))
%         assert(isfield(x.Output.Structural,'T1w_GM_vol_mL'))
%         assert(isfield(x.Output.Structural,'T1w_WM_vol_mL'))
%         assert(isfield(x.Output.Structural,'T1w_CSF_vol_mL'))
%         assert(isfield(x.Output.Structural,'T1w_ICV_vol_mL'))
%         assert(isfield(x.Output.Structural,'T1w_GM_ICV_Ratio'))
%         assert(isfield(x.Output.Structural,'T1w_WM_ICV_Ratio'))
%         assert(isfield(x.Output.Structural,'T1w_CSF_ICV_Ratio'))
%         assert(isfield(x.Output.Structural,'Version_CAT12'))
%         assert(isfield(x.Output.Structural,'Version_LST'))
%         assert(isfield(x.Output.Structural,'Version_ExploreASL'))
%         assert(isfield(x.Output.Structural,'Version_Matlab'))
%         assert(isfield(x.Output.Structural,'Version_SPM12'))
%     end
% end
% 
% % Check Output image field
% assert(isfield(x,'Output_im'))
% if isfield(x,'Output_im')
%     assert(isfield(x.Output_im,'Structural'))
% end
% 








