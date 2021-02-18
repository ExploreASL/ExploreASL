%% BIDS testing script
% Fully test the Flavors by DICOM->BIDS->Legacy import with dedicated results validation
% And then processing all through the ExploreASL pipeline

% Please first run your local path initialization and clone the Flavors Database, the proceed with step
% by step testing
%% Preparation for Henk
clc
pathExploreASL = '/Users/henk/ExploreASL/ExploreASL';
pathTest = '/Users/henk/ExploreASL/ASL/TestBIDS';
cmdCloneFlavors = 'git clone https://github.com/ExploreASL/FlavorDatabase.git';

%% Preparation for Jan
clc
pathExploreASL = '/home/janpetr/code/ExploreASL';
pathTest = '/pet/projekte/asl/data/BIDS';
cmdCloneFlavors = 'git clone git@github.com:ExploreASL/FlavorDatabase.git';

%% Clone the flavors database if necessary
if ~exist(fullfile(pathTest,'FlavorDatabase'), 'dir')
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
%xASL_test_BIDSFlavorsFull(pathExploreASL,pathTest,[0 0 0 0 0 1 0]);

% Check the pipeline results
%xASL_test_BIDSFlavorsFull(pathExploreASL,pathTest,[0 0 0 0 0 0 1]);

% TODO
% - check legacy conversion
% - run the pipeline
% - check the pipeline results
% - convert the datapar.json if provided