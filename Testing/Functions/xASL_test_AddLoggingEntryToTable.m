function logTable = xASL_test_AddLoggingEntryToTable(nameFlavor,logTable,logStruct)
%xASL_test_AddLoggingEntryToTable testing script
%
% FORMAT:       logTable = xASL_test_AddLoggingEntryToTable(nameFlavor,logTable,logStruct)
%
% INPUT:        nameFlavor - Name of the flavor (STRING, REQUIRED)
%               logTable   - Logging table (TABLE, REQUIRED)
%               logStruct  - Logging struct (STRUCT, REQUIRED)
%
% OUTPUT:       logTable   - Logging table
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Add entry to log table.
%
% EXAMPLE:      logTable = xASL_test_AddLoggingEntryToTable(nameFlavor,logTable,logStruct);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL

    % Get number of log entries
    numLogEntries = size(logStruct,2);
    
    % Iterate over log entries
    for iEntry=1:numLogEntries
        thisStruct = logStruct(iEntry);
        thisStruct.name = cellstr(nameFlavor);
        thisStruct.stack = xASL_test_StackToString(thisStruct.stack);
        thisStruct.message = cellstr(thisStruct.message);
        thisStruct.stack = cellstr(thisStruct.stack);
        thisRow = struct2table(thisStruct);
        logTable = [logTable;thisRow];
    end

end