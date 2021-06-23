function xASL_test_Flavors(pathExploreASL, pathTest, bTest, x)
%xASL_test_Flavors Runs the complete testing of Flavors including import from DICOM to BIDS, processing and comparison
%
% FORMAT: xASL_test_Flavors(pathExploreASL, pathTest[, bTest, x])
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
%   x              - x structure (OPTIONAL, DEFAULT = run Initialization)
%
% OUTPUT: 
%  n/a
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Runs the full testing on import and processing of the FlavorsDatabase. The testing directory
%              path has to be provided with the FlavorsDatabase subdirectory containig the Flavors - this 
%              subdirectory is read, but not modified. New directories are created for that inside the test
%              directory.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     xASL_test_Flavors('/home/user/ExploreASL', '/home/user/tmp', [0 1 1 0 1], x)
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

% Change directory to ExploreASL root folder
cd(pathExploreASL);

if nargin < 4 || isempty(x)
    % Remove existing paths
    thisDirectory = pwd;
    if ~isempty(pathExploreASL)
        xASL_adm_RemoveDirectories(pathExploreASL);
        cd(pathExploreASL);
    end
    % Initialize ExploreASL
    ExploreASL_Initialize;
    cd(thisDirectory);
end

% Initialize the working paths
flavorsPath  = fullfile(pathTest, 'FlavorDatabase');
conversionPath = fullfile(pathTest, 'TmpConversion');
referencePath = fullfile(pathTest, 'TmpReference');

% Default dataPar.json for the testing that is fast to run
defaultDataPar.x.subject_regexp = '^sub-.*$';
defaultDataPar.x.bUseMNIasDummyStructural = 1; % when no structural data, use ASL-MNI registration
defaultDataPar.x.Quality = 0; % speed up testing
defaultDataPar.x.DELETETEMP = 1;

%% 1. Make a temporary copy of the Flavors data
if bTest(1)
	% Copy the DICOM source data apart
	fprintf('Copying the source data to a temporary folder: %s\n', conversionPath);
    fprintf('Copying:   ');
	xASL_adm_CreateDir(conversionPath);
    
    if xASL_exist(flavorsPath,'dir')
        fList = xASL_adm_GetFileList(flavorsPath, [], 'List', [], 1);
    else
        error('No FlavorDatabase folder found...');
    end
    fprintf('  '); % Make empty spaces for xASL_TrackProgress
	for iList = 1:length(fList)
        % Display the overall progress and create the conversion directories
		xASL_TrackProgress(iList, length(fList));
		xASL_adm_CreateDir(fullfile(conversionPath,fList{iList}));
		xASL_adm_CreateDir(fullfile(conversionPath,fList{iList}, 'sourcedata'));
		
        % Copy the sourcedata to the temporary conversion folder
		SourceDir = fullfile(flavorsPath, fList{iList}, 'sourcedata');
        DestinationDir = fullfile(conversionPath, fList{iList}, 'sourcedata');
        xASL_Copy(SourceDir, DestinationDir, 1); 
		
        % Copy the JSON meta files to the temporary conversion folder
		xASL_Copy(fullfile(flavorsPath,fList{iList},'sourceStructure.json'), fullfile(conversionPath,fList{iList},'sourceStructure.json'), 1);
		xASL_Copy(fullfile(flavorsPath,fList{iList},'studyPar.json'), fullfile(conversionPath,fList{iList},'studyPar.json'), 1);
		if exist(fullfile(flavorsPath,fList{iList},'dataPar.json'),'file')
			xASL_Copy(fullfile(flavorsPath,fList{iList},'dataPar.json'), fullfile(conversionPath,fList{iList},'dataPar.json'), 1);
		end
	end
	fprintf('\n');
	
	% Create a copy of the reference BIDS data
	fprintf('Copying the BIDS reference data to a temporary folder: %s\n', referencePath);
	xASL_adm_CreateDir(referencePath);
	
	fList = xASL_adm_GetFileList(flavorsPath, [], 'List', [], 1);
    fprintf('  '); % Make empty spaces for xASL_TrackProgress
	for iList = 1:length(fList)
		xASL_TrackProgress(iList,length(fList));
		xASL_adm_CreateDir(fullfile(referencePath,fList{iList}));
		xASL_adm_CreateDir(fullfile(referencePath,fList{iList}, 'rawdata'));

        % Copy rawdata folder for comparison
        xASL_Copy(fullfile(flavorsPath,fList{iList}, 'rawdata'), fullfile(referencePath,fList{iList}, 'rawdata'), 1); 
	end
end

%% 2. Run the conversion of source data to BIDS
if bTest(2)
	xASL_test_Flavors_DCM2BIDS(conversionPath, x);
end

%% 3. Run the comparison of converted BIDS with the reference data
if bTest(3)
	% List all studies in the import directory
	filenameCompare = xASL_adm_GetFileList(conversionPath, '^.+$', 'List', [], true);
	for iCompare = 1:length(filenameCompare)
		% Compare the imported data in the 'rawdata' subdirectory with the counterpart
		fprintf('%s\n', ['Dataset: '  filenameCompare{iCompare}]);
		xASL_bids_CompareStructures(fullfile(conversionPath, filenameCompare{iCompare}, 'rawdata'),...
            fullfile(referencePath, filenameCompare{iCompare}, 'rawdata'),[],[],0,1);
	end
end

%% 4. Run the BIDS to Legacy conversion
if bTest(4)
	% Go through all studies
	ListFolders = xASL_adm_GetFileList(conversionPath, '^.+$', 'FPListRec', [0 Inf], 1);
	for iList=1:numel(ListFolders)
		% Convert only those containing raw data
		if exist(fullfile(ListFolders{iList},'rawdata'),'dir')
			
			% Currently, we clean the old data for unix only
			if isunix
				if exist(fullfile(ListFolders{iList}, 'derivatives'), 'dir')
					diary('off');
					fclose('all'); % ensure that no file is locked
					system(['rm -rf ' fullfile(ListFolders{iList}, 'derivatives')]);
				end
			else
				% Use xASL_delete on windows
                diary('off');
                fclose('all'); % ensure that no file is locked
                xASL_delete(fullfile(ListFolders{iList}, 'derivatives'),true);
			end
			
			% Run the legacy conversion
			% Check if a dataPar is provided, otherwise use the defaults
			fListDataPar = xASL_adm_GetFileList(ListFolders{iList},'(?i)(^dataPar.*\.json$)', 'FPList', [], 0);
			if length(fListDataPar) < 1
				% Fill the dataPars with default parameters
				pathDefaultDataPar = fullfile(ListFolders{iList},'dataPar.json');
				spm_jsonwrite(pathDefaultDataPar,defaultDataPar);
				ExploreASL(ListFolders{iList}, [0 0 0 1], 0, 0);
				xASL_delete(pathDefaultDataPar);
			else
				% Fill the dataPars with the provided parameters
				ExploreASL(ListFolders{iList}, [0 0 0 1], 0, 0);
			end
			
		end
	end
end

%% 5. Run the comparison of data converted to the legacy format with the reference data
if bTest(5)
	error('Not yet implemented');
end

%% 6. Run ExploreASL on all Legacy-converted data
if bTest(6)
	ListFolders = xASL_adm_GetFileList(conversionPath, '^.+$', 'FPListRec', [0 Inf], 1);
	for iList=1:numel(ListFolders)
		pathDerivatives = fullfile(ListFolders{iList},'derivatives');
		% Process data that were converted to derivatives
		if exist(pathDerivatives,'dir')
			pathDerivatives = fullfile(pathDerivatives,'ExploreASL');
			if exist(pathDerivatives,'dir')
				ExploreASL(ListFolders{iList}, 0, [1 1 0], 0); % Don't run population module
			end
		end
	end
end

%% 7. Run the comparison of processed legacy-format data with the reference data
if bTest(7)
	error('Not yet implemented');
end

end