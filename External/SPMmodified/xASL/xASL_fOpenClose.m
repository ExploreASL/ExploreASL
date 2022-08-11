function [fileID, fCloseOutput] = xASL_fOpenClose(pathFile, bOpenClose, bVerbose, bOverwrite)
%xASL_fOpenClose Checks all open files in Matlab and opens or closes only the requested file
%
% FORMAT: [fcloseOutput] = xASL_fOpenClose(pathFile[, bVerbose])
%
% INPUT:
%   pathFile         - path to file (CHAR, REQUIRED)
%   bOpenClose       - true for open and false for close (BOOLEAN, REQUIRED)
%   bVerbose         - Provide verbose output to screen (BOOLEAN, OPTIONAL, DEFAULT = 0)
%   bOverwrite       - When opening a new file, deletes the old file (if
%                      exists) (OPTIONAL, DEFAULT = 1)
%
% OUTPUT:
%   fileID           - fileID (number) for the opened file (DEFAULT = NaN)
%   fcloseOutput     - Report NaN for no closures, 0 for succesful closure, non-zero for errors (OPTIONAL)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:       Close a specified file, if it was open for writing in Matlab.
%                    When Matlab wants to write to a file, it can lock it
%                    for other programs; and when Matlab crashes it can
%                    stay locked. fclose('all') is a universal solution but
%                    closes all open files, which may not be desired when
%                    e.g., logging is done in the background.
% 
%                    Set bOpenClose to open for opening a file if it was not open, or
%                    to close for closing it if it is open
%
% EXAMPLE:
% xASL_fOpenClose('/datasetRoot/ExploreASL/log/importSummary.tsv', 1);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright (c) 2022- ExploreASL


% Admin
if nargin<4 || isempty(bOverwrite)
    bOverwrite = true;
end

if nargin<3 || isempty(bVerbose)
    bVerbose = false;
end
if nargin<2 || isempty(bOpenClose)
    error('Choose open or close');
end
if nargin<1
    error('Input path missing');
elseif isempty(pathFile)
    error('Input path was empty');
elseif exist(pathFile, 'dir')
    error('Input path was a directory but should be a file instead');
end

% defaults
fcloseOutput = NaN;
fileID = NaN;
fileWasOpen = false;

fileIDs = fopen('all'); % Obtain an overview of open files
for iID=1:length(fileIDs)
    fileName{iID} = fopen(fileIDs(iID));
    if strcmp(fileName{iID}, pathFile) % check if the open file is the one requested
        fileWasOpen = true;
        fileID = fileIDs(iID);

        % Always close the file (and then either it's opened again for
        % appending below or not)
        fcloseOutput = fclose(fileIDs(iID)); % close the file for writing

        if fcloseOutput~=0
            warning(['Something went wrong closing ' fileName{iID} ' for writing in Matlab']);

        elseif bVerbose
            fprintf('%s\n', ['Closed ' fileName{iID} ' for writing in Matlab']);
        end
    end
end

% Opening a file, if requested
if bOpenClose
%     if ~fileWasOpen
%         if bOverwrite % make sure to create a new file
%             xASL_delete(pathFile, 0);
%         end
%     end

    fileID = fopen(pathFile, 'at'); % 'r' = reading 'w' writing discarding existing contents
    % a = writing appending to existing contents t = text
end


end




%% OLD
% fileIDs = fopen('all'); % Obtain an overview of open files
% for iID=1:length(fileIDs)
%     fileName{iID} = fopen(fileIDs(iID));
%     if strcmp(fileName{iID}, pathFile) % check if the open file is the one requested
%         fileWasOpen = true;
%         fileID = fileIDs(iID);
% 
%         % Closing a file
%         if ~bOpenClose
%             fcloseOutput = fclose(fileIDs(iID)); % close the file for writing
%     
%             if fcloseOutput~=0
%                 warning(['Something went wrong closing ' fileName{iID} ' for writing in Matlab']);
%     
%             elseif bVerbose
%                 fprintf('%s\n', ['Closed ' fileName{iID} ' for writing in Matlab']);
%             end
%         end
%     end
% end
% 
% % Opening a file
% if ~fileWasOpen
%     if bOverwrite
%         xASL_delete(pathFile, 0);
%     end
% 
%     fileID = fopen(pathFile, 'at'); % 'r' = reading 'w' writing discarding existing contents
%     % a = writing appending to existing contents t = text
% end