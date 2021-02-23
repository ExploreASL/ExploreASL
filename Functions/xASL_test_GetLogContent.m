function [logContent] = xASL_test_GetLogContent(rootDir, printContent, storeRelativePath, exportTable)
%xASL_test_GetLogContent Get warnings and errors from log files
%
% FORMAT: [logContent] = xASL_test_GetLogContent(rootDir, [printContent], [storeRelativePath], [exportTable])
%
% INPUT:
%        rootDir            - Case root directory (REQUIRED, DEFAULT = user input)
%        printContent       - Print warnings and errors (OPTIONAL, DEFAULT = false, BOOLEAN)
%        storeRelativePath  - Store relative path in logContent table instead of module name (OPTIONAL, DEFAULT = true, BOOLEAN)
%        exportTable        - Export a tsv or xlsx file containing the log content (OPTIONAL, DEFAULT = 1, BOOLEAN)
%                             (0 = no export, 1 = TSV export, 2 = XLSX export, 3 = both TSV & XLSX export)
%
% OUTPUT:
%        logContent         - table containing warnings and errors
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      Get warnings and errors from log files.
%
% 0. Input check
% 1. Load all log files
% 2. Iterate over log files
% 3. Optional: Print log content
% 4. Optional: Export (0 = no export, 1 = TSV export, 2 = XLSX export)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          To extract all warnings and errors from all log files
%                   in a directory, you can run this script with the following command.
%                   These settings will not print the warning or error messages, a relative
%                   path will be stored for each file and both TSV and XLSX report files will be generated.
%                   Check out the logContent table if you want to see the results from within Matlab.
%                   rootDir = '.\Test_Runs\TestDataSet';
%                   [logContent] = xASL_test_GetLogContent(rootDir,0,1,3);
%
% REFERENCES:       ...
% __________________________________
% Copyright @ 2015-2021 ExploreASL

    %% -----------------------------------------------------------------------------------------------------------------------------------------------------
    %% 0. Input Check
    if (nargin < 1) || isempty(rootDir)
        if ~usejava('desktop') || ~usejava('jvm') || ~feature('ShowFigureWindows')
            rootDir = input('Insert the root directory: ');
        else
            rootDir = uigetdir(pwd,'Select the root directory...');
        end
    end
    
    if (nargin < 2) || isempty(printContent)
        printContent = false;
    end
    
    if (nargin < 3) || isempty(storeRelativePath)
        storeRelativePath = true;
    end
    
    if (nargin < 4) || isempty(exportTable)
        exportTable = 1;
    end
    
    % Initialize table
    logContent = array2table(zeros(0,7));
    logContent.Properties.VariableNames = {'Module','Subject','Type','Message','File','Line','Content'};
    
    %% -----------------------------------------------------------------------------------------------------------------------------------------------------
    %% 1. Get all log files within root directory
    
    % All log files
    fileList = xASL_adm_GetFileList(rootDir, '^.+\.log$', 'FPListRec');
    
    % Import logs
    fileListImportLogs = xASL_adm_GetFileList(rootDir, '^import_log_.+\.txt$', 'FPListRec');
    
    % Concat lists
    fileList = vertcat(fileList, fileListImportLogs);
    
    % Switch on to check where the script fails
    debugMode = false;
    
    %% -----------------------------------------------------------------------------------------------------------------------------------------------------
    %% 2. Iterate over log files
    for iFile=1:numel(fileList)
        % Get current file path
        curFile = fileList{iFile,1};
        [~, printName] = fileparts(curFile);
        % Print log file
        if debugMode
            fprintf('Log file: %s\n',printName);
        end
        % Fallback subject definition
        curSubject = 'unknown';
        % Current file array
        currentFileArray = strsplit(curFile,filesep);
        % Get subject if there is one
        for level=1:numel(currentFileArray)
            if strfind(lower(currentFileArray{level}),'sub')
                curSubject = currentFileArray{level};
            elseif strfind(lower(currentFileArray{level}),'population')
                curSubject = 'Population';
            end
        end
        
        %% Extract warnings and errors from current log file
        warningsInFile = extractWarnings(curFile,'Warning:','ExploreASL_Master','In <a');
        warningsInFile = getLastFileWarning(warningsInFile,'in <');
        errorsInFile = extractWarnings(curFile,'ERROR: Job iteration terminated!','CONT: but continue with next iteration!');
        errorsInFile = getLastFileError(errorsInFile,'error using <','error in <');
        relativeFileName = strrep(curFile,rootDir,'');
        
        %% Add current warnings
        if ~isempty(warningsInFile{1,1})
            for thisWarning=1:size(warningsInFile,1)
                currentWarning = warningsInFile(thisWarning,1);
                lastWarningFile = warningsInFile(thisWarning,2);
                correspondingLine = warningsInFile(thisWarning,3);
                mainWarnings = warningsInFile(thisWarning,4);
                if storeRelativePath
                    logContent = [logContent;{relativeFileName,curSubject,'Warning',mainWarnings,lastWarningFile,correspondingLine,currentWarning}];
                else
                    logContent = [logContent;{printName,curSubject,'Warning',mainWarnings,lastWarningFile,correspondingLine,currentWarning}];
                end
            end
        end

        %% Add current errors
        if ~isempty(errorsInFile{1,1})
            for thisError=1:size(errorsInFile,1)
                currentError = errorsInFile(thisError,1);
                lastErrorFile = errorsInFile(thisError,2);
                correspondingLine = errorsInFile(thisError,3);
                mainErrors = errorsInFile(thisError,4);
                if storeRelativePath
                    logContent = [logContent;{relativeFileName,curSubject,'Error',mainErrors,lastErrorFile,correspondingLine,currentError}];
                else
                    logContent = [logContent;{printName,curSubject,'Error',mainErrors,lastErrorFile,correspondingLine,currentError}];
                end
            end
        end
    end
    
    %% -----------------------------------------------------------------------------------------------------------------------------------------------------
    %% 3. Optional: Print log content
    if printContent
        fprintf('====================================================================================================\n')
        for curPrint=1:numel(logContent.Content)
            fprintf('Log file: %s\n',logContent.Module{curPrint});
            for linePrint=1:size(logContent.Content{curPrint},1)
                fprintf('%s\n',logContent.Content{curPrint,1}{linePrint});
            end
            fprintf('====================================================================================================\n')
        end
    end
    
    %% -----------------------------------------------------------------------------------------------------------------------------------------------------
    %% 4. Optional: Export (0 = no export, 1 = TSV export, 2 = XLSX export)
    exportLogContent(rootDir,logContent,exportTable);
    
    if isempty(logContent)
        fprintf('No warnings or errors found...\n');
    end
   
end

%% Export files
% INPUT:    rootDir (directory where the resulting files will be stored, CHAR, REQUIRED)
%           logContent (table containing warnings & errors, TABLE, REQUIRED)
%           exportTable (0 = no export, 1 = export TSV, 2 = export XLSX, REQUIRED)
% OUTPUT:   n/a (only stores files)
function exportLogContent(rootDir,logContent,exportTable)

    % Convert warnings & errors from cell to char array
    logContent = logContentCellToChar(logContent);
    if exportTable==1 && ~isempty(logContent)
        % Convert to cell array
        logContentCell = table2cell(logContent);
        % Export TSV file
        try
            if xASL_exist(fullfile(rootDir,'logContent.tsv'))
                fprintf('Deleting previous logContent.tsv file...\n');
                xASL_delete(fullfile(rootDir,'logContent.tsv'));
            end
            xASL_tsvWrite(logContentCell, fullfile(rootDir,'logContent.tsv'), true);
            fprintf('logContent.tsv exported...\n')
        catch ME
            fprintf('%s\n', ME.message);
        end
        if xASL_exist(fullfile(rootDir,'logContent.xlsx'))
            xASL_delete(fullfile(rootDir,'logContent.xlsx'));
            fprintf('logContent.xlsx deleted...\n')
        end
    end
    if exportTable==2 && ~isempty(logContent)
        % Export table
        try
            if xASL_exist(fullfile(rootDir,'logContent.xlsx'))
                fprintf('Deleting previous logContent.xlsx file...\n');
                xASL_delete(fullfile(rootDir,'logContent.xlsx'));
            end
            writetable(logContent,fullfile(rootDir,'logContent.xlsx'),'WriteVariableNames',true);
            fprintf('logContent.xlsx exported...\n')
            if xASL_exist(fullfile(rootDir,'logContent.tsv'))
                xASL_delete(fullfile(rootDir,'logContent.tsv'));
                fprintf('logContent.tsv deleted...\n')
            end
        catch ME
            fprintf('%s\n', ME.message);
        end
    end
    if exportTable==3 && ~isempty(logContent)
        % Convert to cell array
        logContentCell = table2cell(logContent);
        % Export tsv and table
        try
            if xASL_exist(fullfile(rootDir,'logContent.tsv'))
                fprintf('Deleting previous logContent.tsv file...\n');
                xASL_delete(fullfile(rootDir,'logContent.tsv'));
            end
            xASL_tsvWrite(logContentCell, fullfile(rootDir,'logContent.tsv'), true);
            fprintf('logContent.tsv exported...\n')
        catch ME
            fprintf('%s\n', ME.message);
        end
        try
            if xASL_exist(fullfile(rootDir,'logContent.xlsx'))
                fprintf('Deleting previous logContent.xlsx file...\n');
                xASL_delete(fullfile(rootDir,'logContent.xlsx'));
            end
            writetable(logContent,fullfile(rootDir,'logContent.xlsx'),'WriteVariableNames',true);
            fprintf('logContent.xlsx exported...\n')
        catch ME
            fprintf('%s\n', ME.message);
        end
    end
end


%% Get last file of warning or error message (first file in content)
% INPUT:    content (cell array containing a warning or error message, CELL, REQUIRED)
%           identifier (char array containing the identifier of the warning or error, CHAR, REQUIRED)
% OUTPUT:   content (same as input, but re-styled and with additional fields describing the main message etc., CELL)
function content = getLastFileWarning(content,identifier)

    % Iterate over warnings/error messages
    for thisContent=1:numel(content)
        % Get current content
        currentContent = content(thisContent,1);
        content{thisContent,2} = 'unknown';
        content{thisContent,3} = 'unknown';
        content{thisContent,4} = 'unknown';
        found = false;
        % Find and restyle warnings
        for line=1:size(currentContent{1,1},1)
            if line==1
                mainMessage = currentContent{1,1}{line};
                startMainMessage = strfind(mainMessage,'Warning: ');
                if ~isempty(startMainMessage)
                    mainMessage = mainMessage(startMainMessage+length('Warning: '):end);
                    mainMessage = strrep(mainMessage,'[','');
                    mainMessage = strrep(mainMessage,']','');
                    mainMessage = strtrim(mainMessage); % This does not seem to work, yet
                else
                    mainMessage = 'unknown';
                end
                content{thisContent,4} = mainMessage;
                content{thisContent,1}{line} = ''; % Remove this warning text
            else
                content{thisContent,1}{line} = ''; % Remove this warning text
            end
            if ~isempty(strfind(lower(currentContent{1,1}{line}),identifier))
                % Check if HTML style text is in current line
                if ~isempty(strfind(currentContent{1,1}{line},'DocCallback('))
                    % Get Message
                    contentMessage = currentContent{1,1}{line};
                    % Edit style
                    startMessage = strfind(contentMessage,'DocCallback(');
                    contentMessageUpdate = contentMessage(startMessage+length('DocCallback(')+1:end);
                    endMessage = strfind(contentMessageUpdate,',');
                    contentMessageUpdate = contentMessageUpdate(1:endMessage(1)-2);
                    % Line message
                    startLine = strfind(contentMessage,'">line');
                    lineMessageUpdate = contentMessage(startLine+length('">line')+1:end);
                    endLine = strfind(lineMessageUpdate,'</a>');
                    lineMessageUpdate = lineMessageUpdate(1:endLine-1);
                    % Store back
                    content{thisContent,1}{line} = [contentMessageUpdate,', line ',lineMessageUpdate];
                    if ~found
                        content{thisContent,2} = contentMessageUpdate;
                        content{thisContent,3} = lineMessageUpdate;
                        found = true;
                    end
                end
            end
        end
    end
end

%% Get last file of warning or error message (first file in content)
% INPUT:    content (cell array containing a warning or error message, CELL, REQUIRED)
%           identifierA (char array containing the identifier of the warning or error, CHAR, REQUIRED)
%           identifierA (char array containing the identifier of the warning or error, CHAR, REQUIRED)
% OUTPUT:   content (same as input, but re-styled and with additional fields describing the main message etc., CELL)
function content = getLastFileError(content,identifierA,identifierB)

    % Iterate over warnings/error messages
    for thisContent=1:numel(content)
        % Get current content
        currentContent = content(thisContent,1);
        content{thisContent,2} = 'unknown';
        content{thisContent,3} = 'unknown';
        content{thisContent,4} = 'unknown';
        found = false;
        % Find and restyle warnings
        for line=1:size(currentContent{1,1},1)
            if line==1
                mainMessage = currentContent{1,1}{line};
                startMainMessage = strfind(mainMessage,'Warning: ');
                startMainErrorMessage = strfind(mainMessage,'ERROR: ');
                if ~isempty(startMainMessage)
                    mainMessage = mainMessage(startMainMessage+length('Warning: '):end);
                    mainMessage = strrep(mainMessage,'[','');
                    mainMessage = strrep(mainMessage,']','');
                    mainMessage = strtrim(mainMessage); % This does not seem to work, yet
                elseif ~isempty(startMainErrorMessage)
                    mainMessage = mainMessage(startMainErrorMessage+length('ERROR: '):end);
                    mainMessage = strrep(mainMessage,'[','');
                    mainMessage = strrep(mainMessage,']','');
                    mainMessage = strtrim(mainMessage); % This does not seem to work, yet
                else
                    mainMessage = 'unknown';
                end
                content{thisContent,4} = mainMessage;
                content{thisContent,1}{line} = ''; % Remove this warning text
            else
                content{thisContent,1}{line} = ''; % Remove this warning text
            end
            if ~isempty(strfind(lower(currentContent{1,1}{line}),identifierA)) || ~isempty(strfind(lower(currentContent{1,1}{line}),identifierB))
                % Check if HTML style text is in current line
                if ~isempty(strfind(currentContent{1,1}{line},'DocCallback('))
                    % Get Message
                    contentMessage = currentContent{1,1}{line};
                    % Edit style
                    startMessage = strfind(contentMessage,'DocCallback(');
                    contentMessageUpdate = contentMessage(startMessage+length('DocCallback(')+1:end);
                    endMessage = strfind(contentMessageUpdate,',');
                    contentMessageUpdate = contentMessageUpdate(1:endMessage(1)-2);
                    % Line message
                    startLine = strfind(contentMessage,'">line');
                    lineMessageUpdate = contentMessage(startLine+length('">line')+1:end);
                    endLine = strfind(lineMessageUpdate,'</a>');
                    lineMessageUpdate = lineMessageUpdate(1:endLine-1);
                    % Store back
                    content{thisContent,1}{line} = [contentMessageUpdate,', line ',lineMessageUpdate];
                    if ~found
                        content{thisContent,2} = contentMessageUpdate;
                        content{thisContent,3} = lineMessageUpdate;
                        found = true;
                    end
                end
            end
        end
    end
end


%% Convert warnings & errors from cell to char array
% INPUT:    logContent (table containing warnings & errors, TABLE, REQUIRED)
% OUTPUT:   logContent (same as input, but the messages are converted from cell arrays to char arrays, TABLE)
function logContent = logContentCellToChar(logContent)
    % Iterate over rows
    for row=1:size(logContent,1)
        thisContentField = '';
        for thisLine=1:size(logContent.Content{row,1},1)
            if iscell(logContent.Content{row,1})
                if ~isempty(logContent.Content{row,1}{thisLine})
                    thisContentField = [thisContentField,logContent.Content{row,1}{thisLine}];
                    thisContentField = [char(thisContentField),', '];
                end
            else
                % Do not modify custom single line warnings
                thisContentField = logContent.Content{row,1};
            end
        end
        logContent.Content{row,1} = thisContentField;
    end
end

%% Extract warnigns or errors (the alternativeStartIdentifier was necessary for warnings starting without 'Warning:')
function contentInFile = extractWarnings(filePath,startIdentifier,endIdentifier,alternativeStartIdentifier)

    % Check input arguments
    if nargin<4
        alternativeStartIdentifier = [];
    end 

    % Read file
    fileLines = readFileIntoCellArray(filePath);
    
    % Initialize cell array
    contentInFile = cell(1,1);
    id = 1;
    
    % Start and end of current warning or error
    startC = NaN;
    endC = NaN;
    
    % Iterate over lines
    line = 1;
    while line<numel(fileLines)
        
        % Get current line
        curLine = char(fileLines(line,1));
        
        % Check for start of warning or error
        if isempty(alternativeStartIdentifier)
            if contains(curLine, startIdentifier)
                startC = line;

                % Search for end of content
                for subline=startC+1:numel(fileLines)
                    curSubLine = char(fileLines(subline,1));
                    % Check for start of warning or error
                    if contains(curSubLine, endIdentifier)
                        endC = subline;                                     % Current line number
                        currentContentToExtract = fileLines(startC:endC);   % Complete content text
                        startC = NaN;                                       % Reset
                        endC = NaN;                                         % Reset
                        contentInFile(id,1) = {currentContentToExtract};    % Add content to cell array
                        id = id+1;
                        line = subline;                                     % Skip lines
                        break;                                              % Skip this loop
                    end
                end
            end
        else % Handle warnings without 'Warning:' text -> 'In <a'
            if contains(curLine, startIdentifier) || contains(curLine, alternativeStartIdentifier)
                startC = line;

                % Search for end of content
                for subline=startC+1:numel(fileLines)
                    % If the last line was reached, this has to be the endC
                    if subline+1>numel(fileLines)
                        endC = subline;                                     % Current line number
                        currentContentToExtract = fileLines(startC:endC);   % Complete content text
                        startC = NaN;                                       % Reset
                        endC = NaN;                                         % Reset
                        contentInFile(id,1) = {currentContentToExtract};    % Add content to cell array
                        id = id+1;
                        line = subline;                                     % Skip lines
                        break;
                    end
                    % Otherwise ...
                    curSubLine = lower(char(fileLines(subline,1)));         % We have to check this and the following line, because the warning main message can have
                    followingSubLine = lower(char(fileLines(subline+1,1))); % multiple lines and if there is no end identifier, we have to exclude this line as well.
                    % Check for start of warning or error (stop message extraction with end identifier)
                    if contains(curSubLine, lower(endIdentifier)) && ~(contains(followingSubLine, 'in <') || contains(followingSubLine, '> in')) 
                        endC = subline;                                     % Current line number
                        currentContentToExtract = fileLines(startC:endC);   % Complete content text
                        startC = NaN;                                       % Reset
                        endC = NaN;                                         % Reset
                        contentInFile(id,1) = {currentContentToExtract};    % Add content to cell array
                        id = id+1;
                        line = subline;                                     % Skip lines
                        break;                                              % Skip this loop
                    end
                    if strcmp(lower(curLine(1:length('warning: '))),'warning: ')   % Custom single line warnings in import_log
                        endC = subline-1;                                   % Current line number
                        currentContentToExtract = fileLines(startC:endC);   % Complete content text
                        startC = NaN;                                       % Reset
                        endC = NaN;                                         % Reset
                        contentInFile(id,1) = {currentContentToExtract};    % Add content to cell array
                        id = id+1;
                        line = subline-1;                                     % Skip lines
                        break;                                              % Skip this loop
                    end
                    if ~(contains(curSubLine, 'in <') || contains(curSubLine, '> in')) && ...  % stop if 'in <' or '> in' does not occur anymore
                       ~(contains(followingSubLine, 'in <') || contains(followingSubLine, '> in')) 
                        endC = subline-1;                                   % Current line number
                        currentContentToExtract = fileLines(startC:endC);   % Complete content text
                        startC = NaN;                                       % Reset
                        endC = NaN;                                         % Reset
                        contentInFile(id,1) = {currentContentToExtract};    % Add content to cell array
                        id = id+1;
                        line = subline;                                     % Skip lines
                        break;                                              % Skip this loop
                    end
                end
            end
        end
        line = line+1;
    end
end

%% Read text file into cell array
function fileLines = readFileIntoCellArray(filePath)

    % Read file
    fileStr = fileread(filePath);
    fileLines = regexp(fileStr, '\r\n|\r|\n', 'split');
    fileLines = fileLines';

end


