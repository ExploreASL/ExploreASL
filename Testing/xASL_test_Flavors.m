function [flavors, testConfig] = xASL_test_Flavors(testConfig, bOnlyRemoveResults, bRunProcessing)
%xASL_test_Flavors Runs the complete testing of Flavors including import from DICOM to BIDS, processing and comparison
%
% FORMAT: [flavors, testConfig] = xASL_test_Flavors(testConfig, bOnlyRemoveResults, bRunProcessing)
% 
% INPUT:
%   testConfig         - Struct describing the test configuration, or path to a JSON (OPTIONAL, DEFAULT = check this file ExploreASL/Testing/testConfig.json)
%   bOnlyRemoveResults - Set to true if you do not want to run test testing, 
%                        but you want to delete existing test data (BOOLEAN, OPTIONAL, DEFAULT = false) 
%   bRunProcessing     - Run processing (BOOLEAN, OPTIONAL, DEFAULT = false) 
%
% OUTPUT:
%   flavors        - Struct containing the loggingTable and other fields
%   testConfig     - Struct containing all relevant testing fields from the corresponding JSON file
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      Fully test the Flavors by DICOM->BIDS->Legacy import with dedicated
%                   results validation and then processing all through the ExploreASL
%                   pipeline. Please first run your local path initialization and clone the
%                   Flavors Database, the proceed with step by step testing.
%                   1. Clean previous test output
%                   2. Run the DCM->BIDS import
%                   3. Check the DCM->BIDS import results
%                   4. Run BIDS->Legacy import
%                   5. Check the the BIDS->Legacy import results
%                   6. Run ExploreASL on all Legacy-converted data
%                   7. Run the comparison of processed legacy-format data with the reference data
%
% Your testConfig.json could look like this e.g.:
%
% {
%    "pathExploreASL":     "...//ExploreASL",
%    "pathFlavorDatabase": "...//FlavorDatabase",
%    "cmdCloneFlavors":    "git clone https://github.com/ExploreASL/FlavorDatabase.git"
% }
%
% Optionally you can add a list of flavors of interest to your testConfig.json:
%
% {
%    "pathExploreASL":     "...//ExploreASL",
%    "pathFlavorDatabase": "...//FlavorDatabase",
%    "cmdCloneFlavors":    "git clone https://github.com/ExploreASL/FlavorDatabase.git",
%    "flavorList":         ["GE_PCASL_3Dspiral_14.0LX_WIP_1", "Philips_PCASL_2DEPI_3.2.1.1_1", "Siemens_PASL_3DGRASE_E11_1"]
% }
%
% To remove test data from the flavor database you can run:
%
% [flavors, testConfig] = xASL_test_Flavors([],true);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [flavors, testConfig] = xASL_test_Flavors;
%
% __________________________________
% Copyright 2015-2022 ExploreASL


%% 0. Admin and initialization
% Dummy output in case this code crashes
flavors = struct;

% clear screen for easier warning/error tracking
clc;

if nargin<2 || isempty(bOnlyRemoveResults)
	bOnlyRemoveResults = false;
end
if nargin<3 || isempty(bRunProcessing)
	bRunProcessing = false;
end

% Check for testConfig

% Get testing path
pathTesting = fileparts(mfilename('fullpath'));

% Check if testConfig.json needs to be read
if nargin<1 || isempty(testConfig)
	if exist(fullfile(pathTesting,'testConfig.json'),'file')
		testConfig = xASL_io_ReadJson(fullfile(pathTesting,'testConfig.json'));
		if ~(isfield(testConfig,'pathExploreASL') && isfield(testConfig,'pathFlavorDatabase') && isfield(testConfig,'cmdCloneFlavors'))
			fprintf('Please add the correct fields to your testConfig.json...\n');
			return
		end
	else
		fprintf('Please add a testConfig.json to the Testing directory of ExploreASL...\n');
		return
	end
end

if isempty(testConfig.pathExploreASL) || isempty(testConfig.pathFlavorDatabase)
	error('The paths to the code and working directory needs to be specified...');
end

% Clone the flavors database if necessary
cd(testConfig.pathExploreASL);
x = ExploreASL;
if exist(testConfig.pathFlavorDatabase, 'dir')
	cd(testConfig.pathFlavorDatabase);
else
	cd(fileparts(testConfig.pathFlavorDatabase));
	system(testConfig.cmdCloneFlavors);
end

% Load database JSON
if xASL_exist(fullfile(testConfig.pathFlavorDatabase,'flavors.json'),'file')
	databaseInfo = xASL_io_ReadJson(fullfile(testConfig.pathFlavorDatabase,'flavors.json'));
	flavors.data = databaseInfo.flavors;
else
	error('The flavors.json file is missing...');
end
    
% Check for flavor list
if ~isfield(testConfig,'flavorList')
    testConfig.directoryList = xASL_adm_GetFileList(testConfig.pathFlavorDatabase, [], false, [], true);

    % Generate a list of flavors that should be tested from flavors.json
    flavorList = {};
    for iFlavor=1:length(flavors.data)
        if flavors.data(iFlavor).bTesting
            flavorList{end+1,1} = flavors.data(iFlavor).name;
        end
    end
    testConfig.flavorList = flavorList;

    % Check if fields are missing in directory
    testConfig.missingFields = setdiff(testConfig.flavorList, testConfig.directoryList);
    if ~isempty(testConfig.missingFields)
        warning(['Flavors specified in flavors.json are missing in directory'])
        display(testConfig.missingFields);
    end

    % Notify of ignored flavors 
    testConfig.ignoredFields = setdiff(testConfig.directoryList, testConfig.flavorList);
    if ~isempty(testConfig.ignoredFields)
        fprintf('The following flavors will be ignored as specified in flavor.json');
        display(testConfig.ignoredFields);
    end
end    

% Logging table
flavors.loggingTable = array2table(zeros(0,3), 'VariableNames',{'message','stack','name'});
flavors.comparisonTable = array2table(zeros(0,4), 'VariableNames',{'flavor','dataset','name','message'});
    
% Change directory to ExploreASL root folder
cd(testConfig.pathExploreASL);

%% 1. Clean previous test output
% Remove output testdata (derivatives and rawdata), keeps the input and reference (sourcedata, configuration JSONs, rawdataReference, derivativesReference)
xASL_test_Flavors_RemoveExistingTestData(testConfig);

% Stop testing pipeline if we only want to remove test data
if bOnlyRemoveResults
	fclose('all');
	diary off;
	return
end

% Remove empty directories
fprintf('\nChecking for empty folders:   ');
for iList=1:length(testConfig.flavorList)
    xASL_TrackProgress(iList, length(testConfig.flavorList));
    currentFolder = fullfile(testConfig.pathFlavorDatabase, testConfig.flavorList{iList});
    if exist(currentFolder, 'dir') && isempty(xASL_adm_GetFileList(currentFolder, '.*', 'FPListRec'))
        warning(['Deleting empty folder: ' testConfig.flavorList{iList}]);
        xASL_delete(currentFolder, 1);
    end
end


%% 2. Run the DCM->BIDS import
xASL_adm_BreakString('RUN DICOM TO BIDS');
flavors.loggingTable = xASL_test_Flavors_DCM2BIDS(testConfig, x, flavors.loggingTable);

%% 3. Check the DCM->BIDS import results
xASL_adm_BreakString('CHECK THE BIDS CONVERSION');
[flavors,~] = xASL_test_Flavors_Compare(testConfig,flavors,'rawdata','rawdataReference');

%% 4. Run BIDS->Legacy import
xASL_adm_BreakString('RUN BIDS TO LEGACY');
xASL_test_Flavors_BIDS2LEGACY(testConfig);

%% 5. Check the the BIDS->Legacy import results
xASL_adm_BreakString('CHECK THE LEGACY CONVERSION');
[flavors,~] = xASL_test_Flavors_Compare(testConfig,flavors,'derivatives','derivativesReference');

% Already save conversion results and ignore some files before processing
flavors = xASL_test_FlavorsSaveResults(flavors, testConfig);

% Processing
if bRunProcessing
	%% 6. Run ExploreASL on all Legacy-converted data
	xASL_adm_BreakString('RUN PROCESSING PIPELINE');
	flavors.loggingTable = xASL_test_Flavors_ExploreASL(testConfig,flavors.loggingTable);
	
	%% 7. Run the comparison of processed legacy-format data with the reference data
	xASL_adm_BreakString('CHECK OF PROCESSING RESULTS');
	flavors.loggingTable = xASL_test_CheckProcessedFlavors(testConfig,flavors.data,flavors.loggingTable);
	
	% Save all testing results
    flavors = xASL_test_FlavorsSaveResults(flavors, testConfig);
end
 
% Clean-up (file handles etc.)
fclose('all');
diary off;
end

%% Save the test results in a .mat file and ignore log files
function flavors = xASL_test_FlavorsSaveResults(flavors, testConfig)

    % Ignore some files
    flavors = xASL_test_IgnoreFiles(flavors);
    
    % Ignore Acknowledgements & GeneratedBy in dataset_description.json, ASL4D.json, ASL4D_Source.json, M0.json, T1.json, FLAIR.json
    flavors = xASL_test_IgnoreSomeFields(flavors, testConfig);
    
    % Save path
    TimeString = datestr(now,'yyyy-mm-dd_HH_MM');
    TimeString(end-2) = 'h';
    savePathMat = fullfile(testConfig.pathExploreASL,'Testing', [TimeString, '_flavor_results.mat']);
    savePathTSV = fullfile(testConfig.pathExploreASL,'Testing', [TimeString, '_flavor_comparison.tsv']);

	save(savePathMat,'flavors','testConfig');
    if ~isempty(flavors.comparisonTable)
        xASL_tsvWrite(table2cell(flavors.comparisonTable), savePathTSV);
    end 

    % Clear console window
    clc
    
    % Print tables
    fprintf('\n');
    fprintf('[\bCOMPARISON TABLE:]\b\n');
    disp(flavors.comparisonTable);
    fprintf('[\bLOGGING TABLE:]\b\n');
    disp(flavors.loggingTable);
    fprintf('\n');

end

%% Ignore version in dataset_description.json, ASL4D.json, ASL4D_Source.json, M0.json, T1.json, FLAIR.json
function flavors = xASL_test_IgnoreSomeFields(flavors,testConfig)

    % Default
    ignoreRows = [];
    
    % Iterate over table
    for iElement = 1:size(flavors.comparisonTable,1)
        currentElement = flavors.comparisonTable(iElement,:);
        currentFlavor = char(table2cell(currentElement(1,'flavor')));
        currentName = char(table2cell(currentElement(1,'name')));
        currentMessage = char(table2cell(currentElement(1,'message')));
        flavorPath = fullfile(testConfig.pathFlavorDatabase,currentFlavor);
        % Check for different file content in dataset_description files
        if ~isempty(regexpi(currentName,'different file content'))
            % Search for JSON file
            filename = [];
            if ~isempty(regexpi(currentMessage,'dataset_description.json'))
                filename = 'dataset_description.json';
            elseif ~isempty(regexpi(currentMessage,'ASL4D.json'))
                filename = 'ASL4D.json';
            elseif ~isempty(regexpi(currentMessage,'ASL4D_Source.json'))
                filename = 'ASL4D_Source.json';
            elseif ~isempty(regexpi(currentMessage,'M0.json'))
                filename = 'M0.json';
			elseif ~isempty(regexpi(currentMessage,'M0_revPE.json'))
                filename = 'M0_revPE.json';
            elseif ~isempty(regexpi(currentMessage,'T1.json'))
                filename = 'T1.json';
			elseif ~isempty(regexpi(currentMessage,'T2.json'))
				filename = 'T2.json';
            elseif ~isempty(regexpi(currentMessage,'FLAIR.json'))
                filename = 'FLAIR.json';
            end
            % Search for dataset_description.json or other JSON files in derivatives
            ignoreRows = xASL_test_CompareFieldsOfJSON(currentMessage,filename,flavorPath,ignoreRows,iElement);
        end
    end

    % Actually remove the corresponding rows
    flavors.comparisonTable(ignoreRows,:) = [];


end


%% Compare JSON files but ignore the version field
function ignoreRows = xASL_test_CompareFieldsOfJSON(currentMessage,filename,flavorPath,ignoreRows,iElement)

    % We only want to match derivatives JSON files
    if ~isempty(regexpi(currentMessage,'ExploreASL')) && ~isempty(filename)
        startExploreASL = regexpi(currentMessage,'ExploreASL');
        pathA = fullfile(flavorPath,'derivatives',currentMessage(startExploreASL:end));
        pathB = fullfile(flavorPath,'derivativesReference',currentMessage(startExploreASL:end));
        if xASL_exist(pathA,'file') && xASL_exist(pathB,'file')
            % Actual comparison
            jsonA = xASL_io_ReadJson(pathA);
            jsonB = xASL_io_ReadJson(pathB);
            % Get fieldnames
            fieldNamesA = fieldnames(jsonA);
            fieldNamesB = fieldnames(jsonB);
            % Check which fields are shared and which different
            sharedFieldsAB = intersect(fieldNamesB,fieldNamesA);
            % Escape fields
            [jsonA,jsonB] = xASL_test_EscapeToUnix(jsonA,jsonB,sharedFieldsAB);
            % Get differences
            diffSharedFields = xASL_test_CheckSharedJSONFields(jsonA,jsonB,sharedFieldsAB,false);
            % Fields that are in B, but missing in A
            missingFields = setdiff(fieldNamesB,fieldNamesA);
            % We only want to remove a row from the comparison table if there are only difference
            % in the Acknowledgements or GeneratedBy fields of the current JSON. It is important 
            % that there are otherwise no missing files or differences in other fields.
            if isempty(missingFields) && ~diffSharedFields
                ignoreRows = [ignoreRows iElement];
            end
        end
    end

end

%% We need to be sure that there are no differences in the shared fields for our JSON version comparison
function diffSharedFields = xASL_test_CheckSharedJSONFields(jsonA,jsonB,sharedFields,diffSharedFields)

    % Iterate over shared fields
    for iField=1:numel(sharedFields)
        curField = sharedFields{iField};
        % Check type first
        if ~strcmp(class(jsonA.(curField)),class(jsonB.(curField)))
            diffSharedFields = true;
        else
            % Now check the content
            strError = xASL_bids_CompareFieldLists(jsonA,jsonB,sharedFields,{'Acknowledgements','GeneratedBy'});
            if ~isempty(strError)
                diffSharedFields = true;
            end
        end
    end

end

%% Escape to unix
function [jsonA,jsonB] = xASL_test_EscapeToUnix(jsonA,jsonB,sharedFields)

    % Make it windows/unix compatible
    for iField=1:numel(sharedFields)
        curField = sharedFields{iField};
        if ischar(jsonA.(curField))
            jsonA.(curField) = strrep(jsonA.(curField),'\\','/');
            jsonA.(curField) = strrep(jsonA.(curField),'\','/');
        end
        if ischar(jsonB.(curField))
            jsonB.(curField) = strrep(jsonB.(curField),'\\','/');
            jsonA.(curField) = strrep(jsonA.(curField),'\','/');
        end
    end

end
