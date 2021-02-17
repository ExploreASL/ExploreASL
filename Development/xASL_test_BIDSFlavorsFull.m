function xASL_test_BIDSFlavorsFull(pathExploreASL,pathTest,bTest)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs the complete testing of BIDS flavors including import, processing and comparison
%
% FORMAT: xASL_test_BIDSFlavorsFull(pathExploreASL,pathTest)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% 
% INPUT:
%   pathExploreASL - full path to the ExploreASL code (REQUIRED)
%   pathTest       - full path to the testing directory containing the FlavorsDatabase with all the data (REQUIRED)
%   bTest          - an array of booleans specifying which subparts of the test are supposed to be run
%                    1. Make a copy of the flavors data
%                    2. Run the DCM->BIDS import
%                    3. Check the DCM->BIDS import results
%                    4. Run BIDS->Legacy import
%                    5. Check the the BIDS->Legacy import results
%                    6. Run the ExploreASL on all datasets
%                    7. Checks the ExploreASL processing results
%                    (OPTIONAL, DEFAULT = [1 1 1 1 1 1 1])
%
% OUTPUT: 
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Runs the full testing on import and processing of the FlavorsDatabase. The testing directory
%              path has to be provided with the FlavorsDatabase subdirectory containig the Flavors - this 
%              subdirectory is read, but not modified. New directories are created for that inside the test
%              directory.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     xASL_test_BIDSFlavorsFull('/home/user/ExploreASL','/home/user/tmp',[0 1 1 0 1])
% __________________________________
% Copyright 2015-2021 ExploreASL

%% 0. Admin and initialization
if nargin < 2 || isempty(pathExploreASL) || isempty(pathTest)
	error('The path to the code and working directory needs to be specified');
end

if nargin < 3 || isempty(bTest)
	bTest = ones(1,7);
end

if length(bTest) < 7
	bTest(end+1:7) = 1;
end

% Initialize ExploreASL
cd(pathExploreASL);
ExploreASL_Master('',0);

% Initialize the working paths
flavorsPath  = fullfile(pathTest, 'FlavorDatabase');
conversionPath = fullfile(pathTest, 'TmpConversion');
referencePath = fullfile(pathTest, 'TmpReference');

%% 1. Make a temporary copy of the Flavors data
if bTest(1)
	% Copy the DICOM source data apart
	fprintf('Copying the source data to a temporary folder: %s\n',conversionPath);
	if ~exist(conversionPath, 'dir')
		mkdir(conversionPath);
	end
	
	fList = xASL_adm_GetFileList(flavorsPath, [], 'List', [], 1);
	for iList = 1:length(fList)
		if ~exist(fullfile(conversionPath,fList{iList}), 'dir')
			mkdir(fullfile(conversionPath,fList{iList}));
		end
		if ~exist(fullfile(conversionPath,fList{iList}, 'sourcedata'), 'dir')
			mkdir(fullfile(conversionPath,fList{iList}, 'sourcedata'));
		end
		system(['cp -r ' fullfile(flavorsPath,fList{iList},'sourcedata','*') ' ' fullfile(conversionPath,fList{iList},'sourcedata')]);
		system(['cp -r ' fullfile(flavorsPath,fList{iList},'sourceStructure.json') ' ' fullfile(conversionPath,fList{iList},'sourceStructure.json')]);
		system(['cp -r ' fullfile(flavorsPath,fList{iList},'studyPar.json') ' ' fullfile(conversionPath,fList{iList},'studyPar.json')]);
	end
	
	% Create a copy of the reference BIDS data
	fprintf('Copying the BIDS reference data to a temporary folder: %s\n',referencePath);
	if ~exist(referencePath,'dir')
		mkdir(referencePath);
	end
	
	fList = xASL_adm_GetFileList(flavorsPath, [], 'List', [], 1);
	for iList = 1:length(fList)
		if ~exist(fullfile(referencePath,fList{iList}), 'dir')
			mkdir(fullfile(referencePath,fList{iList}));
		end
		if ~exist(fullfile(referencePath,fList{iList}, 'rawdata'), 'dir')
			mkdir(fullfile(referencePath,fList{iList}, 'rawdata'));
		end
		system(['cp -r ' fullfile(flavorsPath,fList{iList}, 'rawdata','*') ' ' fullfile(referencePath,fList{iList},'rawdata')]);
	end
end
%% 2. Run the conversion of source data to BIDS
if bTest(2)
	xASL_bids_TestBidsConversion(conversionPath, referencePath, 1, 0);
end

%% Run the comparison
if bTest(3)
	xASL_bids_TestBidsConversion(conversionPath, referencePath, 0, 1);
end
