function [List] = xASL_adm_CopyMoveFileList(OriDir, DstDir, StrRegExp, bMove, bDir, bRecursive, bOverwrite, bVerbose)
% Searches for files/folders and copy/move them to a new folder, keeping original folder structure intact
%
% FORMAT: [List] = xASL_adm_CopyMoveFileList(OriDir, DstDir, StrRegExp, bMove[, bDir, bRecursive, bOverwrite, bVerbose])
%
% INPUT:
%   OriDir     - Source folder to search through (REQUIRED)
%   DstDir     - Destination folder to copy or move found files/folders to (REQUIRED)
%   StrRegExp  - regular expression used to define files/folders to search for (REQUIRED)
%   bMove      - true for moving found files/folders, false for copying found files/folders (REQUIRED)
%   bDir       - 1 for searching for folders, 0 for searching for files. 2 for searching for both files & folders (OPTIONAL, DEFAULT = 2)
%   bRecursive - true for recursive searching, false for searching in first folder layer only (OPTIONAL, DEFAULT = true)
%   bOverwrite - When true, overwrite if destination exists
%                When false, do nothing (OPTIONAL, DEFAULT = false)
%   bVerbose   - Verbose mode (OPTIONAL, DEFAULT = false)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Moves a file to a file, a file to a directory, or a directory to a directory. 
%              It keeps the initial extensions, no unzipping or zipping after the move. 
%              But it makes sure that only one of {{.nii}} and {{.nii.gz}} exists in the destination directory.
%              Useful to split a large database.
%
% EXAMPLE: xASL_adm_CopyMoveFileList(/mnt/c/Backup/ASL/StudyName/analysis, /mnt/c/Backup/ASL/StudyName/analysis_2, '.*_2.*, 'true); Move all second time point folders & files to a new location, keeping directory structure intact
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright (C) 2020 ExploreASL
%
% 2019-04-01 Henk Mutsaerts. This file was created as a part of ExploreASL
% See LICENSE for details.

%% Manage input arguments
if nargin<4
    error('Not enough input arguments');
end
if nargin<5 || isempty(bDir)
    bDir = 2; % this function only makes sense when run recursively, otherwise a manual cp or mv is faster
end
if nargin<6 || isempty(bRecursive)
    bRecursive = true; % this function only makes sense when run recursively, otherwise a manual cp or mv is faster
end
if nargin<7 || isempty(bOverwrite)
    bOverwrite = false;
end
if nargin<8 || isempty(bVerbose)
    bVerbose = false;
end

if bRecursive
    FPList = 'FPListRec';
else
    FPList = 'FPList';
end

%% Get the list
if bVerbose; fprintf('\n%s\n', 'Creating list to copy/move'); end
if bDir>0
    List{1} = xASL_adm_GetFileList(OriDir, StrRegExp, FPList, [0 Inf], true); % search for folders
else
    List{1} = xASL_adm_GetFileList(OriDir, StrRegExp, FPList, [0 Inf], false); % search for files
end
TotalLength = length(List{1});
if bDir==2
   List{2} = xASL_adm_GetFileList(OriDir, StrRegExp, FPList, [0 Inf], false); % search for files as well
   TotalLength = length(List{1})+length(List{2});
end
   
if TotalLength==0 && bVerbose
    fprintf('%s\n', 'No files/folders found, skipping...');
end

%% Create the new folder, if needed
xASL_adm_CreateDir(DstDir);

%% Start the copying/moving
for iList=1:length(List)
    for iL=1:length(List{iList})
        
        SrcPath = List{iList}{iL};
        DstPath = fullfile(DstDir, SrcPath(length(OriDir)+2:end));        
        
        if ~bVerbose
            if iList==1
                xASL_TrackProgress(iL, TotalLength);
            else
                xASL_TrackProgress(iL+length(List{1}), TotalLength);
            end
        elseif bMove
            fprintf('%s\n',['Move ' SrcPath ' -> ' DstPath]);
        else
            fprintf('%s\n',['Copy ' SrcPath ' -> ' DstPath]);
        end

        if bMove
            xASL_Move(SrcPath, DstPath, bOverwrite, false);
        else
            xASL_Copy(SrcPath, DstPath, bOverwrite, false);
        end
    end
end

end