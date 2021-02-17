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

%% Legacy testing below - to be moved to the dedicated testing function
% 2. BIDS2legacy
% 3. ExploreASL




ListFolders = xASL_adm_GetFileList(pathTest, '^rawdata$', 'FPListRec', [0 Inf], 1);

for iList=1:numel(ListFolders)
    
    %% 2. BIDS2Legacy
    DerivativesDir = fullfile(fileparts(ListFolders{iList}), 'derivatives');
    if ~isunix
        warning('Here we expect a unix-ish system');
    end
    if exist(DerivativesDir, 'dir')
        diary('off');
        fclose('all'); % ensure that no file is locked
        system(['rm -rf ' DerivativesDir]);
    end

    xASL_bids_BIDS2Legacy(ListFolders{iList});
    
    %% 3. Run ExploreASL
    PathDataPar = fullfile(DerivativesDir, 'ExploreASL', 'DataPar.json');
    ExploreASL_Master(PathDataPar, 1, 1, [], [], [1 2]); % don't run population module
    
end