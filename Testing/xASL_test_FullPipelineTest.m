function [flavors, testConfig, logContent] = xASL_test_FullPipelineTest(testConfig,onlyRemoveResults,runProcessing)
%xASL_test_FullPipelineTest BIDS testing script
%
% FORMAT: [flavors, testConfig, logContent] = xASL_test_FullPipelineTest(testConfig,onlyRemoveResults,runProcessing)
% 
% INPUT:
%   testConfig        - Struct describing the test configuration (OPTIONAL, DEFAULT = check for file)
%   onlyRemoveResults - Set to true if you do not want to run test testing, 
%                       but you want to delete existing test data (BOOLEAN, OPTIONAL) 
%   runProcessing     - Run processing (BOOLEAN, DEFAULT=true)
%
% OUTPUT:
%   flavors        - Struct containing the loggingTable and other fields
%   testConfig     - Struct containing all relevant testing fields from the corresponding JSON file
%   logContent     - Struct containing the logged warnings and errors
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
% [flavors, testConfig, logContent] = xASL_test_FullPipelineTest([],true);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [flavors, testConfig, logContent] = xASL_test_FullPipelineTest;
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

    % Get warnings & errors from log files
    [logContent] = xASL_test_GetLogContent(testConfig.pathFlavorDatabase,0,1,2);
    
    % Save all testing results
    flavors = xASL_test_FlavorsSaveResults(flavors, testConfig, logContent);
    
    % Clean-up (file handles etc.)
    fclose('all');
    diary off;


end


%% Save the test results in a .mat file and ignore log files
function flavors = xASL_test_FlavorsSaveResults(flavors, testConfig, logContent)

    % Ignore some files
    flavors = xALS_test_IgnoreFiles(flavors);
    
    % Ignore version in dataset_description.json, ASL4D.json, ASL4D_Source.json, M0.json, T1.json, FLAIR.json
   flavors = xALS_test_IgnoreVersion(flavors, testConfig);
    
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
    fprintf('[\bCOMPARISON TABLE:]\b\n');
    disp(flavors.comparisonTable);
    fprintf('[\bLOGGING TABLE:]\b\n');
    disp(flavors.loggingTable);
    fprintf('\n');

end


%% Ignore version in dataset_description.json, ASL4D.json, ASL4D_Source.json, M0.json, T1.json, FLAIR.json
function flavors = xALS_test_IgnoreVersion(flavors,testConfig)

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
            elseif ~isempty(regexpi(currentMessage,'T1.json'))
                filename = 'T1.json';
            elseif ~isempty(regexpi(currentMessage,'FLAIR.json'))
                filename = 'FLAIR.json';
            end
            % Search for dataset_description.json or other JSON files in derivatives
            if ~isempty(regexpi(currentMessage,'ExploreASL')) && ~isempty(filename)
                pathA = fullfile(flavorPath,'derivatives','ExploreASL',filename);
                pathB = fullfile(flavorPath,'derivativesReference','ExploreASL',filename);
                if xASL_exist(pathA,'file') && xASL_exist(pathB,'file')
                    % Actual comparison
                    jsonA = spm_jsonread(pathA);
                    jsonB = spm_jsonread(pathB);
                    % Get fieldnames
                    fieldNamesA = fieldnames(jsonA);
                    fieldNamesB = fieldnames(jsonB);
                    % Check which fields are shared and which different
                    sharedFieldsAB = intersect(fieldNamesB,fieldNamesA);
                    % Fields that are in B, but missing in A
                    missingFields = setdiff(fieldNamesB,fieldNamesA);
                    % Check that there are no fields missing and the only difference is the version
                    if isempty(missingFields)
                        if ~strcmp(jsonA.GeneratedBy.Version,jsonB.GeneratedBy.Version)
                            ignoreRows = [ignoreRows iElement];
                        end
                    end
                end
            end
        end
    end

    % Actually remove the corresponding rows
    flavors.comparisonTable(ignoreRows,:) = [];


end




