function loggingTable = xASL_test_CheckProcessedFlavors(testConfig,loggingTable)
%xASL_test_CheckProcessedFlavors Check processed flavors
%
% FORMAT: loggingTable = xASL_test_CheckProcessedFlavors(testConfig,loggingTable)
%
% INPUT:        testConfig   - Struct containing test infos (STRUCT, REQUIRED)
%               loggingTable - Logging table (REQUIRED)
%
% OUTPUT:       loggingTable - Logging table
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Check processed flavors.
%
% EXAMPLE:      loggingTable = xASL_test_CheckProcessedFlavors(testConfig,loggingTable);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright (c) 2015-2021 ExploreASL
    
    % Iterate over flavors
    for iList=1:numel(testConfig.flavorList)
        currentFlavor = fullfile(testConfig.pathFlavorDatabase,testConfig.flavorList{iList});
        pathDerivatives = fullfile(currentFlavor,'derivatives');
        % Process data that were converted to derivatives
        if exist(pathDerivatives,'dir')
            pathExploreASLflavor = fullfile(pathDerivatives,'ExploreASL');
            if exist(pathExploreASLflavor,'dir')
                % Check current flavor
                loggingTable = xASL_test_CheckProcessedFlavor(pathExploreASLflavor,currentFlavor,loggingTable);
            end
        end
    end

end


%% Check the current processed flavor, if there is something missing then we return a logging entry
function loggingTable = xASL_test_CheckProcessedFlavor(pathFlavor,currentFlavor,loggingTable)
    
    % Get all files in this flavor directory
    allFiles = xASL_adm_GetFileList(pathFlavor,'.','FPListRec');
    
    % Check list
    checkList{1,1} = 'dataPar.json';
    checkList{2,1} = 'dataset_description.json';
    checkList{3,1} = ['xASL_Report_' '.*' '.pdf'];
    checkList{4,1} = ['catreport_' '.*' '.pdf'];
    checkList{5,1} = ['QC_collection_' '.*' '.json'];
    checkList{6,1} = ['xASL_module_Import_' '.*' '.log'];
    checkList{7,1} = ['xASL_module_Structural_' '.*' '.log'];
    checkList{8,1} = ['xASL_module_ASL_' '.*' '.log'];
    checkList{9,1} = 'CBF.nii.gz';
    
    % Iterate over check list
    for iCheck=1:numel(checkList)
        % Default
        checkFailed = true;
        % Check individual files
        for iFile=1:numel(allFiles)
            if ~isempty(regexpi(allFiles{iFile},checkList{iCheck,1}))
                checkFailed = false;
            end
        end
        % If the check failed we need to create an entry
        if checkFailed
            fprintf(2,'Processed flavor check: failed...\n');
            logEntry.message = ['Missing file ' checkList{iCheck,1} ' in derivatives...'];
            logEntry.stack(1).file = 'xASL_test_CheckProcessedFlavors';
            logEntry.stack(1).name = 'CheckProcessedFlavors';
            logEntry.stack(1).line = 1;
            [~, logEntry.name] = xASL_fileparts(currentFlavor);
            loggingTable = xASL_test_AddLoggingEntryToTable(logEntry.name,loggingTable,logEntry);
            clear logEntry
        end
    end

end


