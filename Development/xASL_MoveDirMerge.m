function xASL_MoveDirMerge(SrcDir, DstDir, bOverwrite, bVerbose)
%xASL_MoveDirMerge Extension on xASL_Move, allowing to merge folders instead of overwrite only
%
% FORMAT: xASL_MoveDirMerge(SrcDir, DstDir[, bOverwrite, bVerbose])
%
% INPUT:
%   SrcPath    - Source file or folder
%   DstPath    - Destination file or folder
%   bOverwrite - When true, overwrite if destination exists
%                When false, do nothing in that case (DEFAULT)
%   bVerbose   - Verbose mode (DEFAULT TRUE)
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Moves a directory to a folder. If the destination folder doesnt exist, this works the same as xASL_Copy.
%              If the destination folder exists, this function will merge the contents rather than delete the existing folder.
%              The bOverwrite and bVerbose settings work identical to xASL_Move
%
% EXAMPLE: xASL_MoveDirMerge('c:\User\path1', 'c:\User\path2'); Merge contents, no overwriting
%          xASL_MoveDirMerge('c:\User\path1', 'c:\User\path2',true);  Merge contents, overwriting
% __________________________________
% Copyright (C) 2015-2019 ExploreASL

%% =========================================
%% Start with the same input argument checks

    % Checks input arguments
    if nargin<2
        error('Please provide both SrcPath & DstPath');
    elseif strcmp(SrcDir, DstDir)
        warning('SrcPath & DstPath were equal, xASL_Move skipped for:');
        fprintf('%s\n', SrcDir);
        return;
    end
    
    if nargin<3 || isempty(bOverwrite)
        bOverwrite = false;
    end
    
	if nargin<4 || isempty(bVerbose)
		bVerbose = true;
    end

    
%% =========================================
%% Start merging

    if ~exist(DstDir,'dir') % if destination dir doesnt exist, simply move
        xASL_Move(SrcDir, DstDir, bOverwrite, bVerbose);
    else
        % if the folder already exists, move its contents)
        ContentList = xASL_adm_GetFileList(SrcDir, '.*', 'List', [0 Inf], true); % for dirs
        for iC=1:length(ContentList)
            
            xASL_Move(fullfile(SrcDir, ContentList{iC}), fullfile(DstDir, ContentList{iC}), bOverwrite, bVerbose);
        end
        ContentList = xASL_adm_GetFileList(SrcDir, '.*', 'List', [0 Inf], false); % for files
        for iC=1:length(ContentList)
            xASL_Move(fullfile(SrcDir, ContentList{iC}), fullfile(DstDir, ContentList{iC}), bOverwrite, bVerbose);
        end
        if bOverwrite && isempty(xASL_adm_GetFileList(SrcDir, '.*', 'FPListRec', [0 Inf]))
            xASL_delete(SrcDir);
        end
    end



end

