function [flavors, testConfig, logContent] = xASL_test_FullPipelineTest(testConfig,onlyRemoveResults)
%xASL_test_FullPipelineTest BIDS testing script
%
% FORMAT: [flavors, testConfig, logContent] = xASL_test_FullPipelineTest(testConfig)
% 
% INPUT:
%   testConfig        - Struct describing the test configuration (OPTIONAL, DEFAULT = check for file)
%   onlyRemoveResults - Set to true if you do not want to run test testing, 
%                       but you want to delete existing test data (BOOLEAN, OPTIONAL) 
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
    
    % Run the pipeline
    flavors = xASL_test_Flavors(testConfig, [0 0 0 0 0 1 0], x, flavors);

    % Check the pipeline results
    % flavors = xASL_test_Flavors(testConfig, [0 0 0 0 0 0 1], x, flavors);
    
    % Ignore some files
    flavors = xALS_test_IgnoreFiles(flavors);

    % Get warnings & errors from log files
    [logContent] = xASL_test_GetLogContent(testConfig.pathFlavorDatabase,0,1,2);
    
    % Clean-up (file handles etc.)
    fclose('all');
    diary off;


end



