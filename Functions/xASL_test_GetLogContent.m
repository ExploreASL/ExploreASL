function [logContent] = xASL_test_GetLogContent(rootDir, printContent, storeFullPath, exportTable)
%xASL_test_GetLogContent Get warnings and errors from log files
%
% FORMAT: [logContent] = xASL_test_GetLogContent(rootDir, [printContent], [storeFullPath], [exportTable])
%
% INPUT:
%        rootDir            - Case root directory (OPTIONAL, DEFAULT = user input)
%        printContent       - Print warnings and errors (OPTIONAL, DEFAULT = false, BOOLEAN)
%        storeFullPath      - Store full path in logContent table instead of module name (OPTIONAL, DEFAULT = false, BOOLEAN)
%        exportTable        - Export a tsv or xlsx file containing the log content (OPTIONAL, DEFAULT = 1, BOOLEAN)
%                             (0 = no export, 1 = TSV export, 2 = XLSX export)
%
% OUTPUT:
%        logContent         - table containing warnings and errors
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      Get warnings and errors from log files.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          % Run the GetLogContent script, do not print the messages, store the full path and export an XLSX
%                   rootDir = '.\Test_Runs\TestDataSet';
%                   [logContent] = xASL_test_GetLogContent(rootDir,0,1,2);
%
% REFERENCES:       ...
% __________________________________
% Copyright @ 2015-2021 ExploreASL

    %% Input Check
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
    
    if (nargin < 3) || isempty(storeFullPath)
        storeFullPath = false;
    end
    
    if (nargin < 4) || isempty(exportTable)
        exportTable = 1;
    end
    
    % Initialize table
    logContent = array2table(zeros(0,6));
    logContent.Properties.VariableNames = {'Module','Subject','Type','File','Line','Content'};
    
    %% Get all log files within root directory
    fileList = xASL_adm_GetFileList(rootDir, '^.+\.log$', 'FPListRec');
    
    % Get Matlab version
    versionYear = version('-release');
    versionYear = str2num(versionYear(1:4));
    if versionYear < 2017
    	fprintf('You are using Matlab %s, we recommend Matlab 2017 or newer...\n',string(versionYear));
    end
    
    %% Iterate over log files
    for iFile=1:numel(fileList)
        % Get current file path
        curFile = fileList{iFile,1};
        [~, printName] = fileparts(curFile);
        % Print log file
        fprintf('Log file: %s\n',printName);
        % Fallback subject definition
        curSubject = 'unknown';
        % Current file array
        currentFileArray = strsplit(curFile,filesep);
        % Get subject if there is one
        for level=1:numel(currentFileArray)
            if versionYear<2017 % Backwards compatibility
                if strfind(lower(currentFileArray{level}),'sub')
                    curSubject = currentFileArray{level};
                elseif strfind(lower(currentFileArray{level}),'population')
                    curSubject = 'Population';
                end
            else
                if contains(lower(currentFileArray{level}), 'sub')
                    curSubject = currentFileArray{level};
                elseif contains(lower(currentFileArray{level}), 'population')
                    curSubject = 'Population';
                end
            end
        end
        
        % Extract warnings and errors from current log file
        warningsInFile = extractWarnings(curFile,'Warning:','ExploreASL_Master');
        warningsInFile = getLastFileWarning(warningsInFile,'> in ',versionYear,'DocCallback',',');
        errorsInFile = extractWarnings(curFile,'ERROR: Job iteration terminated!','CONT: but continue with next iteration!');
        errorsInFile = getLastFileError(errorsInFile,'error using ','error in ',versionYear,'DocCallback',',');
        
        % Add current warnings and errors
        if ~isempty(warningsInFile{1,1})
            for thisWarning=1:size(warningsInFile,1)
                currentWarning = warningsInFile(thisWarning,1);
                lastWarningFile = warningsInFile(thisWarning,2);
                correspondingLine = warningsInFile(thisWarning,3);
                if storeFullPath
                    logContent = [logContent;{curFile,curSubject,'Warning',lastWarningFile,correspondingLine,currentWarning}];
                else
                    logContent = [logContent;{printName,curSubject,'Warning',lastWarningFile,correspondingLine,currentWarning}];
                end
            end
        end
        if ~isempty(errorsInFile{1,1})
            for thisError=1:size(errorsInFile,1)
                currentError = errorsInFile(thisError,1);
                lastErrorFile = errorsInFile(thisError,2);
                correspondingLine = errorsInFile(thisError,3);
                if storeFullPath
                    logContent = [logContent;{curFile,curSubject,'Error',lastErrorFile,correspondingLine,currentError}];
                else
                    logContent = [logContent;{printName,curSubject,'Error',lastErrorFile,correspondingLine,currentError}];
                end
            end
        end
        
    end
    
    %% Optional: Print log content
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
    
    %% Optional: Export (0 = no export, 1 = TSV export, 2 = XLSX export)
    if exportTable==1
        % Convert warnings & errors from cell to char array
        logContent = logContentCellToChar(logContent);
        % Convert to cell array
        logContentCell = table2cell(logContent);
        % Export TSV file
        xASL_tsvWrite(logContentCell, fullfile(rootDir,'logContent.tsv'), true);
    elseif exportTable==2
        % Convert warnings & errors from cell to char array
        logContent = logContentCellToChar(logContent);
        % Export table
        writetable(logContent,fullfile(rootDir,'logContent.xlsx'),'WriteVariableNames',true);
    end
   
end

% Get last file of warning or error message (first file in content)
function content = getLastFileWarning(content,identifier,versionYear,startMessage,endMessage)

    % Iterate over warnings/error messages
    for thisContent=1:numel(content)
        % Get current content
        currentContent = content(thisContent,1);
        content{thisContent,2} = 'unkown';
        content{thisContent,3} = 'unkown';
        found = false;
        % Find "first" file
        for line=1:size(currentContent{1,1},1)
            if versionYear<2017 % Backwards compatibility
                if ~isempty(strfind(lower(currentContent{1,1}{line}),identifier)) && ~found
                    contentMessage = currentContent{1,1}{line};
                    startFile = strfind(contentMessage,startMessage)+13;
                    contentMessage = contentMessage(startFile:end);
                    endFile = strfind(contentMessage,endMessage);
                    startLineMessage = strfind(contentMessage,'line ');
                    lineInfo = contentMessage(startLineMessage:end-5);
                    contentMessage = contentMessage(1:endFile(1)-2);
                    content{thisContent,2} = contentMessage;
                    content{thisContent,3} = lineInfo;
                    found = true;
                end
            else
                if contains(lower(currentContent{1,1}{line}),identifier) && ~found
                    contentMessage = currentContent{1,1}{line};
                    startFile = strfind(contentMessage,startMessage)+13;
                    contentMessage = contentMessage(startFile:end);
                    endFile = strfind(contentMessage,endMessage);
                    startLineMessage = strfind(contentMessage,'line ');
                    lineInfo = contentMessage(startLineMessage:end-5);
                    contentMessage = contentMessage(1:endFile(1)-2);
                    content{thisContent,2} = contentMessage;
                    content{thisContent,3} = lineInfo;
                    found = true;
                end
                % Find problems
                if line==size(currentContent{1,1},1) && ~found
                    disp('test');
                    
                end
            end
        end
    end
end

% Get last file of warning or error message (first file in content)
function content = getLastFileError(content,identifierA,identifierB,versionYear,startMessage,endMessage)

    % Iterate over warnings/error messages
    for thisContent=1:numel(content)
        % Get current content
        currentContent = content(thisContent,1);
        content{thisContent,2} = 'unkown';
        content{thisContent,3} = 'unkown';
        found = false;
        % Find "first" file
        for line=1:size(currentContent{1,1},1)
            if versionYear<2017 % Backwards compatibility
                if (~isempty(strfind(lower(currentContent{1,1}{line}),identifierA)) || ~isempty(strfind(lower(currentContent{1,1}{line}),identifierB))) && ~found
                    contentMessage = currentContent{1,1}{line};
                    startFile = strfind(contentMessage,startMessage)+13;
                    contentMessage = contentMessage(startFile:end);
                    endFile = strfind(contentMessage,endMessage);
                    startLineMessage = strfind(contentMessage,'line ');
                    lineInfo = contentMessage(startLineMessage:end-5);
                    contentMessage = contentMessage(1:endFile(1)-2);
                    content{thisContent,2} = contentMessage;
                    content{thisContent,3} = lineInfo;
                    found = true;
                end
            else
                if (contains(lower(currentContent{1,1}{line}),identifierA) || contains(lower(currentContent{1,1}{line}),identifierB)) && ~found
                    contentMessage = currentContent{1,1}{line};
                    startFile = strfind(contentMessage,startMessage)+13;
                    contentMessage = contentMessage(startFile:end);
                    endFile = strfind(contentMessage,endMessage);
                    startLineMessage = strfind(contentMessage,'line ');
                    lineInfo = contentMessage(startLineMessage:end-5);
                    contentMessage = contentMessage(1:endFile(1)-2);
                    content{thisContent,2} = contentMessage;
                    content{thisContent,3} = lineInfo;
                    found = true;
                end
            end
        end
    end
end


% Convert warnings & errors from cell to char array
function logContent = logContentCellToChar(logContent)
    % Iterate over rows
    for row=1:size(logContent,1)
        thisContentField = '';
        for thisLine=1:size(logContent.Content{row,1},1)
            thisContentField = strcat(thisContentField,char(logContent.Content{row,1}(thisLine)));
            thisContentField = strcat(thisContentField,' - ');
        end
        logContent.Content{row,1} = thisContentField;
    end
end

% Extract warnigns or errors
function contentInFile = extractWarnings(filePath,startIdentifier,endIdentifier)

    % Read file
    fileLines = readFileIntoCellArray(filePath);
    
    % Initialize cell array
    contentInFile = cell(1,1);
    id = 1;
    
    % Start and end of current warning or error
    startC = NaN;
    endC = NaN;
    
    % Iterate over lines
    for line=1:numel(fileLines)
        
        % Get current line
        curLine = char(fileLines(line,1));
        
        % Check for start of warning or error
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
                    break;                                              % Skip this loop
                end
            end
        end
    end
end

% Read text file into cell array
function fileLines = readFileIntoCellArray(filePath)

    % Read file
    fileStr = fileread(filePath);
    fileLines = regexp(fileStr, '\r\n|\r|\n', 'split');
    fileLines = fileLines';

end


