function [logContent] = xASL_qc_GetLogContent(rootDir, printContent, storeFullPath, exportTSV)
%xASL_qc_GetLogContent Get warnings and errors from log files
%
% FORMAT: [logContent] = xASL_qc_GetLogContent(rootDir, [printContent], [storeFullPath], [exportTSV])
%
% INPUT:
%        rootDir            - Case root directory (OPTIONAL, DEFAULT = user input)
%        printContent       - Print warnings and errors (OPTIONAL, DEFAULT = false)
%        storeFullPath      - Store full path in logContent table instead of module name (OPTIONAL, DEFAULT = false)
%        exportTSV          - Export a tsv file containing the log content (OPTIONAL, DEFAULT = false)
%
% OUTPUT:
%        logContent         - table containing warnings and errors
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      Get warnings and errors from log files.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          rootDir = '.\Test_Runs\TestDataSet';
%                   [logContent] = xASL_qc_GetLogContent(rootDir,false,false,true);
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
    
    if (nargin < 4) || isempty(exportTSV)
        exportTSV = false;
    end
    
    % Initialize table
    logContent = array2table(zeros(0,4));
    logContent.Properties.VariableNames = {'Module','Subject','Type','Content'};
    
    %% Get all log files within root directory
    fileList = dir(fullfile(rootDir, '**\*.log'));  % Get list of log files and folders in any subfolder
    fileList = fileList(~[fileList.isdir]);         % Remove folders from list
    
    %% Iterate over log files
    for iFile=1:numel(fileList)
        % Get current file path
        curFile = fullfile(fileList(iFile).folder,fileList(iFile).name);
        % Print log file
        fprintf('Log file: %s\n',fileList(iFile).name);
        % Fallback subject definition
        curSubject = 'unknown';
        % Current file array
        currentFileArray = strsplit(curFile,filesep);
        % Get subject if there is one
        for level=1:numel(currentFileArray)
            if contains(lower(currentFileArray(level)), 'sub')
                curSubject = currentFileArray(level);
            elseif contains(lower(currentFileArray(level)), 'population')
                curSubject = 'Population';
            end
        end
        
        % Extract warnings and errors from current log file
        warningsInFile = extractWarnings(curFile,'Warning:','ExploreASL_Master');
        errorsInFile = extractWarnings(curFile,'ERROR: Job iteration terminated!','CONT: but continue with next iteration!');
        
        % Add current warnings and errors
        if ~isempty(warningsInFile{1,1})
            for thisWarning=1:numel(warningsInFile)
                currentWarning = warningsInFile(thisWarning,1);
                if storeFullPath
                    logContent = [logContent;{curFile,curSubject,'Warning',currentWarning}];
                else
                    logContent = [logContent;{fileList(iFile).name,curSubject,'Warning',currentWarning}];
                end
            end
        end
        if ~isempty(errorsInFile{1,1})
            for thisError=1:numel(errorsInFile)
                currentError = errorsInFile(thisError,1);
                if storeFullPath
                    logContent = [logContent;{curFile,curSubject,'Error',currentError}];
                else
                    logContent = [logContent;{fileList(iFile).name,curSubject,'Error',currentError}];
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
    
    %% Optional: Export TSV
    if exportTSV
        % Convert warnings & errors from cell to char array
        logContent = logContentCellToChar(logContent);
        % Convert to cell array
        logContentCell = table2cell(logContent);
        % Export TSV file
        xASL_tsvWrite(logContentCell, fullfile(rootDir,'logContent.tsv'), true);
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


