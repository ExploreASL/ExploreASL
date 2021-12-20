function [flavors, testConfig] = xASL_test_FullPipelineTest(testConfig, onlyRemoveResults, runProcessing)
%xASL_test_FullPipelineTest BIDS testing script
%
% FORMAT: [flavors, testConfig] = xASL_test_FullPipelineTest(testConfig, onlyRemoveResults, runProcessing)
% 
% INPUT:
%   testConfig        - Struct describing the test configuration (OPTIONAL, DEFAULT = check for file)
%   onlyRemoveResults - Set to true if you do not want to run test testing, 
%                       but you want to delete existing test data (BOOLEAN, OPTIONAL, DEFAULT = false) 
%   runProcessing     - Run processing (BOOLEAN, OPTIONAL, DEFAULT = true) 
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
%
% Your testConfig.json could look like this e.g.:
%
% {
%    "pathExploreASL":     "...//ExploreASL",
%    "pathFlavorDatabase": "...//FlavorDatabase",
%    "cmdCloneFlavors":    "git clone git@github.com:ExploreASL/FlavorDatabase.git"
% }
%
% Optionally you can add a list of flavors of intereset to your testConfig.json:
%
% {
%    "pathExploreASL":     "...//ExploreASL",
%    "pathFlavorDatabase": "...//FlavorDatabase",
%    "cmdCloneFlavors":    "git clone git@github.com:ExploreASL/FlavorDatabase.git",
%    "flavorList":         ["GE_PCASL_3Dspiral_14.0LX_1", "Philips_PCASL_2DEPI_3.2.1.1_1", "Siemens_PASL_3DGRASE_E11_1"]
% }
%
% To remove test data from the flavor database you can run:
%
% [flavors, testConfig] = xASL_test_FullPipelineTest([],true);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [flavors, testConfig] = xASL_test_FullPipelineTest;
%
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% Initialization
    if nargin<2 || isempty(onlyRemoveResults)
        onlyRemoveResults = false;
    end
    if nargin<3 || isempty(runProcessing)
        runProcessing = true;
    end

    %% Check for testConfig
    
    % Get testing path
    pathTesting = fileparts(mfilename('fullpath'));
    
    % Check if testConfig.json needs to be read
    if nargin<1 || isempty(testConfig)
        if exist(fullfile(pathTesting,'testConfig.json'),'file')
            testConfig = spm_jsonread(fullfile(pathTesting,'testConfig.json'));
            if ~(isfield(testConfig,'pathExploreASL') && isfield(testConfig,'pathFlavorDatabase') && isfield(testConfig,'cmdCloneFlavors'))
                fprintf('Please add the correct fields to your testConfig.json...\n');
                return
            end
        else
            fprintf('Please add a testConfig.json to the Testing directory of ExploreASL...\n');
            return
        end
    end
    
    
    %% Clone the flavors database if necessary
    cd(testConfig.pathExploreASL);
    x = ExploreASL;
    if exist(testConfig.pathFlavorDatabase, 'dir')
        cd(testConfig.pathFlavorDatabase);
    else
        cd(fileparts(testConfig.pathFlavorDatabase));
        system(testConfig.cmdCloneFlavors);
    end
    
    % Check for flavor list
    if ~isfield(testConfig,'flavorList')
        testConfig.flavorList = xASL_adm_GetFileList(testConfig.pathFlavorDatabase, [], false, [], true);
    end
    
    
    %% Logging table
    flavors.loggingTable = array2table(zeros(0,3), 'VariableNames',{'message','stack','name'});
    flavors.comparisonTable = array2table(zeros(0,4), 'VariableNames',{'flavor','dataset','name','message'});
    

    %% Test execution

    % Remove exiting test data
    flavors = xASL_test_Flavors(testConfig, [1 0 0 0 0 0 0], x, flavors);
    
    % Stop testing pipeline if we only want to remove test data
    if onlyRemoveResults
        logContent = struct;
        fclose('all');
        diary off;
        return
    end

    % Convert to BIDS
    flavors = xASL_test_Flavors(testConfig, [0 1 0 0 0 0 0], x, flavors);

    % Check the BIDS conversion
    flavors = xASL_test_Flavors(testConfig, [0 0 1 0 0 0 0], x, flavors);

    % Convert BIDS to Legacy
    flavors = xASL_test_Flavors(testConfig, [0 0 0 1 0 0 0], x, flavors);

    % Check the Legacy conversion
    flavors = xASL_test_Flavors(testConfig, [0 0 0 0 1 0 0], x, flavors);
    
    % Already save conversion results and ignore some files before processing
    flavors = xASL_test_FlavorsSaveResults(flavors, testConfig);
    
    % Processing
    if runProcessing
    
        % Run the pipeline
        flavors = xASL_test_Flavors(testConfig, [0 0 0 0 0 1 0], x, flavors);

        % Check the pipeline results
        flavors = xASL_test_Flavors(testConfig, [0 0 0 0 0 0 1], x, flavors);
        
    end
    
    % Save all testing results
    flavors = xASL_test_FlavorsSaveResults(flavors, testConfig);
    
    % Clean-up (file handles etc.)
    fclose('all');
    diary off;


end


%% Save the test results in a .mat file and ignore log files
function flavors = xASL_test_FlavorsSaveResults(flavors, testConfig, logContent)

    % Ignore some files
    flavors = xALS_test_IgnoreFiles(flavors);
    
    % Ignore Acknowledgements & GeneratedBy in dataset_description.json, ASL4D.json, ASL4D_Source.json, M0.json, T1.json, FLAIR.json
    flavors = xALS_test_IgnoreSomeFields(flavors, testConfig);
    
    % Save path
    savePath = fullfile(testConfig.pathExploreASL,'Testing','results.mat');
    
    % Check if there is a logContent
    if nargin < 3
        save(savePath,'flavors','testConfig');
    else
        save(savePath,'flavors','testConfig','logContent');
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
function flavors = xALS_test_IgnoreSomeFields(flavors,testConfig)

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
            ignoreRows = xALS_test_CompareFieldsOfJSON(currentMessage,filename,flavorPath,ignoreRows,iElement);
        end
    end

    % Actually remove the corresponding rows
    flavors.comparisonTable(ignoreRows,:) = [];


end


%% Compare JSON files but ignore the version field
function ignoreRows = xALS_test_CompareFieldsOfJSON(currentMessage,filename,flavorPath,ignoreRows,iElement)

    % We only want to match derivatives JSON files
    if ~isempty(regexpi(currentMessage,'ExploreASL')) && ~isempty(filename)
        startExploreASL = regexpi(currentMessage,'ExploreASL');
        pathA = fullfile(flavorPath,'derivatives',currentMessage(startExploreASL:end));
        pathB = fullfile(flavorPath,'derivativesReference',currentMessage(startExploreASL:end));
        if xASL_exist(pathA,'file') && xASL_exist(pathB,'file')
            % Actual comparison
            jsonA = spm_jsonread(pathA);
            jsonB = spm_jsonread(pathB);
            % Get fieldnames
            fieldNamesA = fieldnames(jsonA);
            fieldNamesB = fieldnames(jsonB);
            % Check which fields are shared and which different
            sharedFieldsAB = intersect(fieldNamesB,fieldNamesA);
            % Escape fields
            [jsonA,jsonB] = xASL_test_EscapeToUnix(jsonA,jsonB,sharedFieldsAB);
            % Get differences
            diffSharedFields = xALS_test_CheckSharedJSONFields(jsonA,jsonB,sharedFieldsAB,false);
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
function diffSharedFields = xALS_test_CheckSharedJSONFields(jsonA,jsonB,sharedFields,diffSharedFields)

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




