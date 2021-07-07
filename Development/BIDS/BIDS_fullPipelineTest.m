%% BIDS testing script
% Fully test the Flavors by DICOM->BIDS->Legacy import with dedicated results validation
% And then processing all through the ExploreASL pipeline

% Please first run your local path initialization and clone the Flavors Database, the proceed with step
% by step testing

%% Get Username
clc
if isunix
    [~,username] = system('id -u -n');
    username=username(1:end-1);
else
    username = getenv('username');
end

%% Preparation for Henk
if strcmp(username,'henk')
    pathExploreASL = '/Users/henk/ExploreASL/ExploreASL';
    pathTest = '/Users/henk/ExploreASL';
    cmdCloneFlavors = 'git clone https://github.com/ExploreASL/FlavorDatabase.git';
end

%% Preparation for Jan
if strcmp(username,'janpetr')
	if ismac
		pathExploreASL = '/Volumes/janpetr/ExploreASL/ExploreASL';
		pathTest = '/Volumes/pet/projekte/asl/data/BIDS';
	else
		pathExploreASL = '/home/janpetr/ExploreASL/ExploreASL';
		pathTest = '/pet/projekte/asl/data/BIDS';
	end
	cmdCloneFlavors = 'git clone git@github.com:ExploreASL/FlavorDatabase.git';
end


%% Preparation for Michael
if strcmp(username,'matlab')
	pathExploreASL = 'M:\SoftwareDevelopment\MATLAB\m.stritt\Server_xASL\ExploreASL';
	pathTest = 'M:\SoftwareDevelopment\MATLAB\m.stritt\Server_xASL\FlavorTests';
	cmdCloneFlavors = 'git clone git@github.com:ExploreASL/FlavorDatabase.git';
end

%% Preparation for Beatriz
if strcmp(username, 'bestevespadrela')
	pathExploreASL = '/home/bestevespadrela/ExploreASL/';
	pathTest = '/home/bestevespadrela/lood_storage/divi/Projects/ExploreASL/TestBIDS';
	cmdCloneFlavors = 'git clone git@github.com:ExploreASL/FlavorDatabase.git';
end

%% Preparation for Mathijs
if strcmp(username,'mathijs')
	pathExploreASL = '/home/mdijsselhof/ExploreASL/';
	pathTest = '/home/mdijsselhof/Test_BIDS';
	cmdCloneFlavors = 'git clone git@github.com:ExploreASL/FlavorDatabase.git';
end

%% Clone the flavors database if necessary
cd(pathExploreASL);
[x] = ExploreASL_Initialize; % Initialize ExploreASL
cd(pathTest);
if ~exist(fullfile(pathTest,'FlavorDatabase'), 'dir')
	xASL_adm_CreateDir(pathTest);
	system(cmdCloneFlavors);
end

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
