function [logContent] = xASL_qc_GetLogContent(rootDir)
%xASL_qc_GetLogContent Get warnings and errors from log files
%
% FORMAT: [logContent] = xASL_qc_GetLogContent(rootDir)
%
% INPUT:
%        rootDir            - Case root directory
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
%                   [logContent] = xASL_qc_GetLogContent(rootDir);
%
% REFERENCES:       ...
% __________________________________
% Copyright @ 2015-2021 ExploreASL

    %% Input Check
    if nargin < 1
        error('Missing xASL directory...');
    end
    
    % Initialize table
    logContent = array2table(zeros(0,3));
    logContent.Properties.VariableNames = {'Module','Type','Content'};
    
    %% Get all log files within root directory
    fileList = dir(fullfile(rootDir, '**\*.log'));  % Get list of log files and folders in any subfolder
    fileList = fileList(~[fileList.isdir]);         % Remove folders from list
    
    %% Iterate over log files
    for iFile=1:numel(fileList)
        % Get current file path
        curFile = fullfile(fileList(iFile).folder,fileList(iFile).name);
        % Print log file
        fprintf('Log file: %s\n',fileList(iFile).name);
        
        % Extract warnings and errors from current log file
        warningsInFile = extractWarnings(curFile,'Warning:','ExploreASL_Master');
        errorsInFile = extractWarnings(curFile,'ERROR: Job iteration terminated!','CONT: but continue with next iteration!');
        
        % Add current warnings and errors
        if ~isempty(warningsInFile{1,1})
            for thisWarning=1:numel(warningsInFile)
                currentWarning = warningsInFile(thisWarning,1);
                logContent = [logContent;{fileList(iFile).name,'Warning',currentWarning}];
            end
        end
        if ~isempty(errorsInFile{1,1})
            for thisError=1:numel(errorsInFile)
                currentWarning = errorsInFile(thisError,1);
                logContent = [logContent;{fileList(iFile).name,'Error',currentWarning}];
            end
        end
        
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


