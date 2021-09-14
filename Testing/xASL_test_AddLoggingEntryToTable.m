%% Add entry to log table
function logTable = xASL_test_AddLoggingEntryToTable(nameFlavor,logTable,logStruct)

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