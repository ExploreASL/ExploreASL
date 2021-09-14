function [flavors, testConfig, logContent] = xASL_test_FullPipelineTest
%xASL_test_FullPipelineTest BIDS testing script
%
% FORMAT: [flavors, testConfig, logContent] = xASL_test_FullPipelineTest
% 
% INPUT:
%   n/a
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
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [flavors, testConfig, logContent] = xASL_test_FullPipelineTest;
%
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% Check for testConfig
    pathTesting = fileparts(mfilename('fullpath'));
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

    
    %% Clone the flavors database if necessary
    cd(testConfig.pathExploreASL);
    x = ExploreASL;
    if exist(testConfig.pathFlavorDatabase, 'dir')
        cd(testConfig.pathFlavorDatabase);
    else
        cd(fileparts(testConfig.pathFlavorDatabase));
        system(testConfig.cmdCloneFlavors);
    end
    
    
    %% Logging table
    flavors.loggingTable = array2table(zeros(0,3), 'VariableNames',{'message','stack','name'});
    

    %% Test execution

    % Remove exiting test data
    flavors = xASL_test_Flavors(testConfig.pathExploreASL, testConfig.pathFlavorDatabase, [1 0 0 0 0 0 0], x, flavors);

    % Convert to BIDS
    flavors = xASL_test_Flavors(testConfig.pathExploreASL, testConfig.pathFlavorDatabase, [0 1 0 0 0 0 0], x, flavors);

    % Check the BIDS conversion
    flavors = xASL_test_Flavors(testConfig.pathExploreASL, testConfig.pathFlavorDatabase, [0 0 1 0 0 0 0], x, flavors);

    % Convert BIDS to Legacy
    flavors = xASL_test_Flavors(testConfig.pathExploreASL, testConfig.pathFlavorDatabase, [0 0 0 1 0 0 0], x, flavors);

    % Check the Legacy conversion
    flavors = xASL_test_Flavors(testConfig.pathExploreASL, testConfig.pathFlavorDatabase, [0 0 0 0 1 0 0], x, flavors);

    % Run the pipeline
    flavors = xASL_test_Flavors(testConfig.pathExploreASL, testConfig.pathFlavorDatabase, [0 0 0 0 0 1 0], x, flavors);

    % Check the pipeline results
    % flavors = xASL_test_Flavors(testConfig.pathExploreASL, testConfig.pathFlavorDatabase, [0 0 0 0 0 0 1], x, flavors);

    % Get warnings & errors from log files
    [logContent] = xASL_test_GetLogContent(testConfig.pathFlavorDatabase,0,1,2);
    
    % Clean-up (file handles etc.)
    fclose('all');
    diary off;


end



