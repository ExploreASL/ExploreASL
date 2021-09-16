function flavors = xASL_test_Flavors(testConfig, bTest, x, flavors)
%xASL_test_Flavors Runs the complete testing of Flavors including import from DICOM to BIDS, processing and comparison
%
% FORMAT: flavors = xASL_test_Flavors(pathExploreASL, pathFlavorDatabase[, bTest, x],flavors)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% 
% INPUT:
%   testConfig         - struct which contains the paths to the ExploreASL
%                        installation and the testing/flavor repository (REQUIRED) 
%   bTest              - an array of booleans specifying which subparts of the test are supposed to be run
%                        1. Make a copy of the flavors data
%                        2. Run the DCM->BIDS import
%                        3. Check the DCM->BIDS import results
%                        4. Run BIDS->Legacy import
%                        5. Check the the BIDS->Legacy import results
%                        6. Run the ExploreASL on all datasets
%                        7. Checks the ExploreASL processing results
%                        (OPTIONAL, DEFAULT = [1 1 1 1 1 1 1])
%   x                  - x structure (OPTIONAL, DEFAULT = run Initialization)
%   flavors            - struct containing flavor related fields (RECOMMENDED, STRUCT)
%
% OUTPUT: 
%   flavors            - struct containing flavor related fields
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Runs the full testing on import and processing of the FlavorsDatabase. The testing directory
%              path has to be provided with the FlavorsDatabase subdirectory containig the Flavors - this 
%              subdirectory is read, but not modified. New directories are created for that inside the test
%              directory.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     flavors = xASL_test_Flavors(testConfig, [1 0 0 0 0 0 0], x, flavors);
%
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% 0. Admin and initialization
    if isempty(testConfig.pathExploreASL) || isempty(testConfig.pathFlavorDatabase)
        error('The paths to the code and working directory needs to be specified...');
    end
    if nargin < 3 || isempty(bTest)
        bTest = ones(1,7);
    end
    if length(bTest) < 7
        bTest(end+1:7) = 1;
    end
    if nargin < 4
        % The flavors struct is normally being defined in xASL_test_FullPipelineTest
        warning('We recommend to provide the flavors struct based on xASL_test_FullPipelineTest...');
        flavors.loggingTable = array2table(zeros(0,3), 'VariableNames',{'message','stack','name'});
        flavors.comparisonTable = array2table(zeros(0,4), 'VariableNames',{'flavor','dataset','name','message'});
    end

    % Change directory to ExploreASL root folder
    cd(testConfig.pathExploreASL);

    if nargin < 4 || isempty(x)
        % Remove existing paths
        thisDirectory = pwd;
        if ~isempty(testConfig.pathExploreASL)
            xASL_adm_RemoveDirectories(testConfig.pathExploreASL);
            cd(testConfig.pathExploreASL);
        end
        % Initialize ExploreASL
        ExploreASL_Initialize;
        cd(thisDirectory);
    end
    
    % Check for flavor list
    if ~isfield(testConfig,'flavorList')
        testConfig.flavorList = xASL_adm_GetFileList(testConfig.pathFlavorDatabase, [], false, [], true);
    end
    

    %% 1. Remove existing test data
    if bTest(1)
        xASL_test_Flavors_RemoveExistingTestData(testConfig);
    end
    

    %% 2. Run the conversion of source data to BIDS
    if bTest(2)
        flavors.loggingTable = xASL_test_Flavors_DCM2BIDS(testConfig, x, flavors.loggingTable);
    end
    

    %% 3. Run the comparison of converted BIDS with the reference data
    if bTest(3)
        fprintf('\n================================== CHECK THE BIDS CONVERSION =================================\n');
        [flavors,~] = xASL_test_Flavors_Compare(testConfig,flavors,'rawdata','rawdataReference');
    end
    

    %% 4. Run the BIDS to Legacy conversion
    if bTest(4)
        xASL_test_Flavors_BIDS2LEGACY(testConfig);
    end
    

    %% 5. Run the comparison of data converted to the legacy format with the reference data
    if bTest(5)
        fprintf('\n================================= CHECK THE LEGACY CONVERSION ================================\n');
        [flavors,~] = xASL_test_Flavors_Compare(testConfig,flavors,'derivatives','derivativesReference');
    end
    

    %% 6. Run ExploreASL on all Legacy-converted data
    if bTest(6)
        flavors.loggingTable = xASL_test_Flavors_ExploreASL(testConfig,flavors.loggingTable);
    end
    

    %% 7. Run the comparison of processed legacy-format data with the reference data
    if bTest(7)
        error('Not yet implemented...');
    end
    

end



