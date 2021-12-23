function loggingTable = xASL_test_CheckProcessedFlavors(testConfig, flavorData, loggingTable)
%xASL_test_CheckProcessedFlavors Check processed flavors
%
% FORMAT: loggingTable = xASL_test_CheckProcessedFlavors(testConfig, flavorData, loggingTable)
%
% INPUT:        testConfig   - Struct containing test infos (STRUCT, REQUIRED)
%               flavorData   - Reference data of the flavors (STRUCT ARRAY, REQUIRED)
%               loggingTable - Logging table (REQUIRED)
%
% OUTPUT:       loggingTable - Logging table
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Check processed flavors.
%
% EXAMPLE:      loggingTable = xASL_test_CheckProcessedFlavors(testConfig, flavorData, loggingTable);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright (c) 2015-2021 ExploreASL
    
    % Iterate over flavors
    for iList=1:numel(testConfig.flavorList)
        currentFlavor = fullfile(testConfig.pathFlavorDatabase,testConfig.flavorList{iList});
        currentData = xASL_test_GetCurrentFlavorData(flavorData,currentFlavor);
        pathDerivatives = fullfile(currentFlavor,'derivatives');
        % Process data that were converted to derivatives
        if exist(pathDerivatives,'dir')
            pathExploreASLflavor = fullfile(pathDerivatives,'ExploreASL');
            if exist(pathExploreASLflavor,'dir')
                if ~isempty(currentData)
                    % Check current flavor
                    loggingTable = xASL_test_CheckProcessedFlavor(pathExploreASLflavor,currentFlavor,loggingTable,currentData);
                else
                    warning('No reference data for this flavor...');
                end
            else
                warning('No ExploreASL data...');
            end
        else
            warning('No derivatives data...');
        end
    end

end


%% Check the current processed flavor, if there is something missing then we return a logging entry
function loggingTable = xASL_test_CheckProcessedFlavor(pathFlavor,currentFlavor,loggingTable,currentData)
    
    % Get all files in this flavor directory
    allFiles = xASL_adm_GetFileList(pathFlavor,'.','FPListRec');
    
    % Check list
    it = 1;
    checkList{it,1} = 'dataPar.json'; it=it+1;
    checkList{it,1} = 'dataset_description.json'; it=it+1;
    checkList{it,1} = ['xASL_Report_' '.*' '.pdf']; it=it+1;
    if ~currentData.dummyStructural
        checkList{it,1} = ['catreport_' '.*' '.pdf']; it=it+1;
    end
    checkList{it,1} = ['QC_collection_' '.*' '.json']; it=it+1;
    checkList{it,1} = ['xASL_module_Import_' '.*' '.log']; it=it+1;
    checkList{it,1} = ['xASL_module_Structural_' '.*' '.log']; it=it+1;
    checkList{it,1} = ['xASL_module_ASL_' '.*' '.log']; it=it+1;
    checkList{it,1} = 'CBF.nii.gz'; it=it+1;
    
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


%% Get the data of this flavor
function currentData = xASL_test_GetCurrentFlavorData(flavorData,currentFlavor)

    currentData = [];
    for iFlavor = 1:numel(flavorData)
        [~, currentName] = xASL_fileparts(currentFlavor);
        if strcmp(flavorData(iFlavor).name,currentName)
            currentData = flavorData(iFlavor);
            return 
        end
    end

end


