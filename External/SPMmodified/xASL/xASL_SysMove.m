function xASL_SysMove(SrcPath, DstPath, bForce)
%xASL_SysMove Move a file or a directory. Manage different OSes
%
% FORMAT: xASL_SysMove(SrcPath, DstPath[, bForce])
%
% INPUT:
%   SrcPath    - Source file or folder
%   DstPath    - Destination file or folder
%   bForce     - When true, overwrite if destination exists
%                When false, do nothing in that case (DEFAULT)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Moves a file to a file, a file to a directory, or a directory to a directory. SBypass inefficient matlab stuff on linux and windows, but can only move on same file system!
% EXAMPLE: xASL_SysMove('c:\User\path\file.nii', 'c:\User\path2\file.nii'); No overwriting
%          xASL_SysMove('c:\User\path\file.nii', 'c:\User\path2\file.nii',true);  Overwriting
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright (C) 2015-2019 ExploreASL

    %% Manage input arguments
    if nargin<2
        error('Please provide both SrcPath & DstPath');
    elseif strcmp(SrcPath,DstPath)
        warning('SrcPath & DstPath were equal, no action required');
        fprintf('%s\n', SrcPath);
        return;
    end
    if nargin<3 || isempty(bForce)
        bForce = false;
    end
	if ~xASL_exist(SrcPath,'file') % Check if source exists
		error(['xASL_Move: Source doesn''t exist: ' SrcPath]);
    end

    %% Start moving
    if isunix
        SrcPath = xASL_adm_UnixPath(SrcPath);
        DstPath = xASL_adm_UnixPath(DstPath);
        
        if bForce
            strforce = '-f ';
        else
            strforce = [];
        end
        [status,result] = system(['mv ' strforce SrcPath ' ' DstPath]);
        if status~=0
            error('xASL_Move: Error moving %s to %s: %s', SrcPath, DstPath, result);
        end
    elseif ispc
        if bForce
            strforce = '/Y';
        else
            strforce = [];
        end
        [status,result] = system(['move ' strforce ' "' SrcPath '" "' DstPath '"']);
        if status~=0
            error('xASL_Move: Error moving %s to %s: %s', SrcPath, DstPath, result);
        end
    else
        if bForce
            movefile(SrcPath, DstPath, 'f');
        else
            movefile(SrcPath, DstPath);
        end
    end
    
    
end