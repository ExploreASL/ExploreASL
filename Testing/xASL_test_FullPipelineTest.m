%% xASL_test_FullPipelineTest: BIDS testing script
%
% Fully test the Flavors by DICOM->BIDS->Legacy import with dedicated
% results validation and then processing all through the ExploreASL
% pipeline. Please first run your local path initialization and clone the
% Flavors Database, the proceed with step by step testing.
%
% Your testConfig.json could look like this e.g.:
%
% {
%    "pathExploreASL": "...//ExploreASL",
%    "pathTest": "...//FlavorTests",
%    "cmdCloneFlavors": "git clone git@github.com:ExploreASL/FlavorDatabase.git"
% }


%% Check for testConfig
pathTesting = fileparts(mfilename('fullpath'));
pathExploreASL = fileparts(pathTesting);

if exist(fullfile(pathTesting,'testConfig.json'),'file')
    testConfig = spm_jsonread(fullfile(pathTesting,'testConfig.json'));
    if isfield(testConfig,'pathExploreASL') && isfield(testConfig,'pathTest') && isfield(testConfig,'cmdCloneFlavors')
        pathExploreASL = testConfig.pathExploreASL;
        pathTest = testConfig.pathTest;
        cmdCloneFlavors = testConfig.cmdCloneFlavors;
    else
        fprintf('Please add the correct fields to your testConfig.json...\n');
        return
    end
else
    fprintf('Please add a testConfig.json to the Testing directory of ExploreASL...\n');
    return
end

%% Clone the flavors database if necessary
cd(pathExploreASL);
x = ExploreASL;
cd(pathTest);
if ~exist(fullfile(pathTest,'FlavorDatabase'), 'dir')
	xASL_adm_CreateDir(pathTest);
	system(cmdCloneFlavors);
end

%% Remove old testing data
% If these folders contain anything more than only the designated tests,
% this confuses the comparisons, so these have to be deleted first.

% Initialize the working paths
flavorsPath  = fullfile(pathTest, 'FlavorDatabase');
conversionPath = fullfile(pathTest, 'TmpConversion');
referencePath = fullfile(pathTest, 'TmpReference');

xASL_delete(conversionPath, true);
xASL_delete(referencePath, true);

%% Test execution

% Prepare the data
xASL_test_Flavors(pathExploreASL, pathTest, [1 0 0 0 0 0 0], x);

% Convert to BIDS
xASL_test_Flavors(pathExploreASL, pathTest, [0 1 0 0 0 0 0], x);

% Check the BIDS conversion
xASL_test_Flavors(pathExploreASL, pathTest, [0 0 1 0 0 0 0], x);

% Convert BIDS to Legacy
xASL_test_Flavors(pathExploreASL, pathTest, [0 0 0 1 0 0 0], x);

% Check the Legacy conversion
%xASL_test_Flavors(pathExploreASL, pathTest, [0 0 0 0 1 0 0], x);

% Run the pipeline
xASL_test_Flavors(pathExploreASL, pathTest, [0 0 0 0 0 1 0], x);

% Check the pipeline results
%xASL_test_Flavors(pathExploreASL, pathTest, [0 0 0 0 0 0 1], x);

% Get warnings & errors from log files
[logContent] = xASL_test_GetLogContent(pathTest,0,1,2);

% TODO
% - check legacy conversion
% - check the pipeline results

