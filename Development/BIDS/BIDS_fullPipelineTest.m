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
    pathTest = '/Users/henk/ExploreASL/ASL/TestBIDS';
    cmdCloneFlavors = 'git clone https://github.com/ExploreASL/FlavorDatabase.git';
end

%% Preparation for Jan
if strcmp(username,'janpetr')
    pathExploreASL = '/home/janpetr/ExploreASL/ExploreASL';
    pathTest = '/pet/projekte/asl/data/BIDS';
    cmdCloneFlavors = 'git clone git@github.com:ExploreASL/FlavorDatabase.git';
end

%% Preparation for Michael
if strcmp(username,'matlab')
    pathExploreASL = 'M:\SoftwareDevelopment\MATLAB\m.stritt\ExploreASL';
    pathTest = 'M:\SoftwareDevelopment\MATLAB\m.stritt\TestBIDS';
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
if ~exist(fullfile(pathTest,'FlavorDatabase'), 'dir')
    cd(pathExploreASL);
    ExploreASL_Master('',0); % initialize ExploreASL
    xASL_adm_CreateDir(pathTest);
    cd(pathTest);
    system(cmdCloneFlavors);
end

%% Test execution

% Prepare the data
xASL_test_BIDSFlavorsFull(pathExploreASL,pathTest,[1 0 0 0 0 0 0]);

% Convert to BIDS
xASL_test_BIDSFlavorsFull(pathExploreASL,pathTest,[0 1 0 0 0 0 0]);

% Check the BIDS conversion
xASL_test_BIDSFlavorsFull(pathExploreASL,pathTest,[0 0 1 0 0 0 0]);

% Convert BIDS to Legacy
xASL_test_BIDSFlavorsFull(pathExploreASL,pathTest,[0 0 0 1 0 0 0]);

% Check the Legacy conversion
%xASL_test_BIDSFlavorsFull(pathExploreASL,pathTest,[0 0 0 0 1 0 0]);

% Run the pipeline
xASL_test_BIDSFlavorsFull(pathExploreASL,pathTest,[0 0 0 0 0 1 0]);

% Check the pipeline results
%xASL_test_BIDSFlavorsFull(pathExploreASL,pathTest,[0 0 0 0 0 0 1]);

[logContent] = xASL_test_GetLogContent(pathTest,0,1,2);

% TODO
% - check legacy conversion
% - check the pipeline results

