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
                logEntry = xASL_test_GetProcessedEntry(pathExploreASLflavor);
                % Add the log entry if it does exist
                loggingTable = xASL_test_AddLoggingEntry(loggingTable,logEntry);
            end
        end
    end

end


%% Check the current processed flavor, if there is something missing then we return a logging entry
function logEntry = xASL_test_GetProcessedEntry(pathFlavor)

    % Default empty entry
    logEntry = [];

    % ...

end


%% Add the logging entry
function loggingTable = xASL_test_AddLoggingEntry(loggingTable,logEntry)

    if ~isempty(logEntry)

    end

end


