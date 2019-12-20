function xASL_SysCopy(SrcPath, DstPath, bOverwrite, bVerbose)
%xASL_SysCopy Copies a file or a directory. Manage different OSes
%
% FORMAT: xASL_SysCopy(SrcPath, DstPath, bOverwrite, bVerbose)
%
% INPUT:
%   SrcPath    - Source file or folder
%   DstPath    - Destination file or folder
%   bForce     - When true, overwrite if destination exists
%                When false, do nothing in that case (DEFAULT)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Copies a file to a file or a directory to a directory. Bypass inefficient matlab stuff on linux and windows,
%              but can only move on the same file system.
% EXAMPLE: xASL_SysCopy('c:\User\path\file.nii', 'c:\User\path2\file.nii'); No overwriting
%          xASL_SysCopy('c:\User\path\file.nii', 'c:\User\path2\file.nii',true);  Overwriting
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright (C) 2015-2019 ExploreASL

    %% Check input parameters
	if nargin<3 || isempty(bOverwrite)
		bOverwrite = false;
    end
	if nargin<4 || isempty(bVerbose)
		bVerbose = false;
    end
    if strcmp(SrcPath,DstPath)
        return; % skip copying
    end

    %% Do the copying

    if exist(SrcPath, 'dir') && exist(DstPath, 'dir') && bOverwrite
        if bVerbose
            warning('xASL_Copy:We merge folders now, be careful, this could go wrong when e.g. NIfTIs exist both unzipped & zipped');
        end
    end

    if isunix()
        if bOverwrite
            unix(['cp -r -f ' SrcPath ' ' DstPath]); % -n is short for --noclober
        else
            unix(['cp -r ' SrcPath ' ' DstPath]);
        end
    else
        if bOverwrite
            copyfile(SrcPath, DstPath, 'f');
        else
            copyfile(SrcPath, DstPath);
        end
    end
end


