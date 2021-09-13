function xASL_test_Flavors(pathExploreASL, pathFlavorDatabase, bTest, x)
%xASL_test_Flavors Runs the complete testing of Flavors including import from DICOM to BIDS, processing and comparison
%
% FORMAT: xASL_test_Flavors(pathExploreASL, pathFlavorDatabase[, bTest, x])
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% 
% INPUT:
%   pathExploreASL     - full path to the ExploreASL code (REQUIRED)
%   pathFlavorDatabase - full path to the FlavorDatabase directory with all the data (REQUIRED)
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
    if nargin < 2 || isempty(pathExploreASL) || isempty(pathFlavorDatabase)
        error('The paths to the code and working directory needs to be specified...');
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


    %% Default dataPar.json for the testing that is fast to run
    
    % Default regular expression
    defaultDataPar.x.dataset.subjectRegexp = '^sub-.*$';
    
    % When there is no structural data, use ASL-MNI registration
    defaultDataPar.x.modules.asl.bUseMNIasDummyStructural = 1;
    
    % Speed up testing and delete temporary files
    defaultDataPar.x.settings.Quality = 0; 
    defaultDataPar.x.settings.DELETETEMP = 1;


    %% 1. Remove existing test data
    if bTest(1)
        xASL_test_Flavors_RemoveExistingTestData(pathFlavorDatabase);
    end
    

    %% 2. Run the conversion of source data to BIDS
    if bTest(2)
        xASL_test_Flavors_DCM2BIDS(pathFlavorDatabase, x);
    end
    

    %% 3. Run the comparison of converted BIDS with the reference data
    if bTest(3)
        flavors.listNII2BIDS = xASL_test_Flavors_Compare_NII2BIDS(pathFlavorDatabase);
    end
    

    %% 4. Run the BIDS to Legacy conversion
    if bTest(4)
        xASL_test_Flavors_BIDS2LEGACY(pathFlavorDatabase,defaultDataPar);
    end
    

    %% 5. Run the comparison of data converted to the legacy format with the reference data
    if bTest(5)
        flavors.listBIDS2LEGACY = xASL_test_Flavors_Compare_BIDS2LEGACY(pathFlavorDatabase);
    end
    

    %% 6. Run ExploreASL on all Legacy-converted data
    if bTest(6)
        flavors.loggingTable = xASL_test_Flavors_ExploreASL(pathFlavorDatabase);
    end
    

    %% 7. Run the comparison of processed legacy-format data with the reference data
    if bTest(7)
        error('Not yet implemented...');
    end
    

end


%% Remove existing test data
function xASL_test_Flavors_RemoveExistingTestData(pathFlavorDatabase)

    fprintf('Remove existing test data...\n');

    % Load the list of the directories
    flavorList = xASL_adm_GetFileList(pathFlavorDatabase, [], false, [], true);
    for iFlavor = 1:length(flavorList)
        % Get current flavor directory
        thisFlavorDir = fullfile(pathFlavorDatabase,flavorList{iFlavor});

        % Check for temp folder
        if xASL_exist(fullfile(thisFlavorDir,'temp'), 'dir')==7
            xASL_delete(fullfile(thisFlavorDir,'temp'),1);
        end
        
        % Check for rawdata folder
        if xASL_exist(fullfile(thisFlavorDir,'rawdata'), 'dir')==7
            xASL_delete(fullfile(thisFlavorDir,'rawdata'),1);
        end

        % Check for derivatives folder
        if xASL_exist(fullfile(thisFlavorDir,'derivatives'), 'dir')==7
            xASL_delete(fullfile(thisFlavorDir,'derivatives'),1);
        end

        % Check for log and other files
        allFilesInThisDir = xASL_adm_GetFileList(thisFlavorDir,'^.+$');

        % Do not delete sourceStructure.json, studyPar.json and dataPar.json
        for iFile=1:numel(allFilesInThisDir)
            [~, thisFile, thisExtension] = fileparts(allFilesInThisDir{iFile});
            if ~strcmp([thisFile thisExtension],'sourceStructure.json') && ...
                    ~strcmp([thisFile thisExtension],'studyPar.json') && ...
                    ~strcmp([thisFile thisExtension],'dataPar.json')
                xASL_delete(allFilesInThisDir{iFile},1);
            end
        end
    end


end


%% Compare NII2BIDS conversion
function flavorList = xASL_test_Flavors_Compare_NII2BIDS(pathFlavorDatabase)

    % List all studies in the import directory
    flavorList = xASL_adm_GetFileList(pathFlavorDatabase, '^.+$', 'List', [], true);
    
    % Iterate over flavors
    for iCompare = 1:length(flavorList)
        
        % Compare the imported data in the 'rawdata' subdirectory with the counterpart
        fprintf('%s\n', ['Dataset: '  flavorList{iCompare}]);
        convertedRawData = fullfile(pathFlavorDatabase, flavorList{iCompare}, 'rawdata');
        referenceRawData = fullfile(pathFlavorDatabase, flavorList{iCompare}, 'rawdataReference');
        
        % Run comparison script
        [flavorList{iCompare,2},flavorList{iCompare,3}] = ...
            xASL_bids_CompareStructures(referenceRawData,convertedRawData,[],[],0,1);
    end

end


%% Run the BIDS to Legacy conversion
function xASL_test_Flavors_BIDS2LEGACY(pathFlavorDatabase,defaultDataPar)

    % Go through all studies
    ListFolders = xASL_adm_GetFileList(pathFlavorDatabase, '^.+$', 'FPListRec', [0 Inf], 1);
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


%% Run ExploreASL on all Legacy-converted data
function loggingTable = xASL_test_Flavors_ExploreASL(pathFlavorDatabase)

    % Logging table
    loggingTable = array2table(zeros(0,2), 'VariableNames',{'message','stack'});

    % Get flavor list
    ListFolders = xASL_adm_GetFileList(pathFlavorDatabase, '^.+$', 'FPListRec', [0 Inf], 1);
    
    % Iterate over flavors
    for iList=1:numel(ListFolders)
        pathDerivatives = fullfile(ListFolders{iList},'derivatives');
        % Process data that were converted to derivatives
        if exist(pathDerivatives,'dir')
            pathDerivatives = fullfile(pathDerivatives,'ExploreASL');
            if exist(pathDerivatives,'dir')
                % Don't run population module
                xFlavor = ExploreASL(ListFolders{iList}, 0, [1 1 0], 0);
                if isfield(xFlavor,'logging')
                    loggingTable = xASL_test_AddLoggingEntryToTable(loggingTable,xFlavor.logging);
                end
            end
        end
    end

end


%% Run the comparison of data converted to the legacy format with the reference data
function flavorList = xASL_test_Flavors_Compare_BIDS2LEGACY(pathFlavorDatabase)


    % List all studies in the import directory
    flavorList = xASL_adm_GetFileList(pathFlavorDatabase, '^.+$', 'List', [], true);
    
    % Iterate over flavors
    for iCompare = 1:length(flavorList)
        
        % Compare the imported data in the 'rawdata' subdirectory with the counterpart
        fprintf('%s\n', ['Dataset: '  flavorList{iCompare}]);
        convertedDerivativesData = fullfile(pathFlavorDatabase, flavorList{iCompare}, 'derivatives');
        referenceDerivativesData = fullfile(pathFlavorDatabase, flavorList{iCompare}, 'derivativesReference');
        
        % Run comparison script
        [flavorList{iCompare,2},flavorList{iCompare,3}] = ...
            xASL_bids_CompareStructures(referenceDerivativesData,convertedDerivativesData,[],[],0,1);
    end


end


%% Add entry to log table
function logTable = xASL_test_AddLoggingEntryToTable(logTable,logStruct)

    % Get number of log entries
    numLogEntries = size(logStruct,2);
    
    % Iterate over log entries
    for iEntry=1:numLogEntries
        thisStruct = logStruct(iEntry);
        thisStruct.stack = xASL_test_StackToString(thisStruct.stack);
        thisStruct.message = cellstr(thisStruct.message);
        thisStruct.stack = cellstr(thisStruct.stack);
        thisRow = struct2table(thisStruct);
        logTable = [logTable;thisRow];
    end

end






