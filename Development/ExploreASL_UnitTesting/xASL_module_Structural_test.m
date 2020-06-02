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

% Get ExploreASL directory
xASLdir = uigetdir(pwd, 'Select ExploreASL directory...');

% Define test directory
testDir = uigetdir(pwd, 'Select testing directory...');
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
assert(isfield(x,'out'))                 % Check if new ... field was created
% ...


% WORK IN PROGRESS
% Which fields always have to be in the resulting x structure?



% Run again and also save 'result' variable to assert this too
% assert(isfield(result,''))















