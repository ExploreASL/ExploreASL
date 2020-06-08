%xASL_Initialize_test Script to test the xASL_Initialize function
%
% FORMAT:       RESULT = runtests('xASL_Initialize_test');
% 
% INPUT:        None
%
% OUTPUT:       Console window
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This script is supposed to run various unit tests of xASL_Initialize:
%
%           1) Run a test using the default TestDataSet inputs with low quality setting
%           2) Run a test using the default TestDataSet inputs with high quality setting
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLES: RESULT = runtests('xASL_Initialize_test');
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

% Field names
LIST =     {'name', 'subject_regexp', 'Quality', 'bNativeSpaceAnalysis', 'nWorkers', 'iWorker', 'S',...
            'ProcessData', 'MyPath', 'D', 'stopaftererrors', 'dryrun', 'bOverwrite', 'Version', 'Q',...
            'SpaghettiDir', 'HistogramDir', 'StatsMaps', 'SPMDIR', 'SPMpath', 'SPMVERSION', 'skull',...
            'BILAT_FILTER', 'WBmask', 'bReproTesting', 'Pediatric_Template', 'bLesionFilling', 'bAutoACPC',...
            'Segment_SPM12', 'M0_conventionalProcessing', 'bGetControlLabelOrder', 'TotalSubjects',...
            'nTotalSubjects', 'exclusion', 'SESSIONS', 'nSessions', 'TimePointTotalSubjects',...
            'ExcludedSubjects', 'SUBJECTS', 'TotalInclusionList', 'nSubjects', 'nExcluded', 'nSubjectsSessions',...
            'nTimePointsTotal', 'nTimePointTotalSubjects', 'TimePointSubjects', 'nTimePoints',...
            'nTimePointSubjects', 'TimePointExcluded', 'nTimePointExcluded', 'DELETETEMP', 'P'};

% Field types
TYPE_LIST =   { 'char', 'char', 'double', 'double', 'double', 'double', 'struct', 'logical', 'char',...
                'struct', 'double', 'logical', 'logical', 'char', 'struct', 'char', 'char', 'char',...
                'char', 'char', 'char', 'logical', 'logical', 'logical', 'logical', 'logical', 'logical',...
                'logical', 'logical', 'logical', 'logical', 'cell', 'double', 'cell', 'cell', 'double',...
                'cell', 'char', 'cell', 'double', 'double', 'double', 'double', 'double', 'double',...
                'cell', 'double', 'double', 'cell', 'double', 'double', 'struct'};


%% Test 1: Default TestDataSet with low quality setting
fprintf('Initialize test input...\n')
DataParPath = fullfile(testDir,'TestFolder','TestDataSet','DataParameters_LowQ.json');
ProcessData = true;
iWorker = 1;
nWorkers = 1;

% Run test
x = ExploreASL_Initialize(DataParPath, ProcessData, iWorker, nWorkers);

% Check quality setting
assert(x.Quality==0)

% Use assert to check if all fields exist
for i=1:length(LIST)
    assert(isfield(x,LIST{i}))
end

% Check field formats
for i=1:length(LIST)
    % Compare format of current field with list
    assert(strcmp(class(x.(LIST{i})),TYPE_LIST{i}))
end


%% Test 2: Default TestDataSet with high quality setting
fprintf('Initialize test input...\n')
DataParPath = fullfile(testDir,'TestFolder','TestDataSet','DataParameters_HiQ.json');
ProcessData = true;
iWorker = 1;
nWorkers = 1;

% Run test
x = ExploreASL_Initialize(DataParPath, ProcessData, iWorker, nWorkers);

% Check quality setting
assert(x.Quality==1)

% Use assert to check if all fields exist
for i=1:length(LIST)
    assert(isfield(x,LIST{i}))
end

% Check field formats
for i=1:length(LIST)
    % Compare format of current field with list
    assert(strcmp(class(x.(LIST{i})),TYPE_LIST{i}))
end








