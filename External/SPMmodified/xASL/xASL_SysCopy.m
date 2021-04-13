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
        warning('SrcPath & DstPath were equal, no action required');
        fprintf('%s\n', SrcPath);
        return; % skip copying
    end
    
    
    %% Do the copying

    if exist(SrcPath, 'dir') && exist(DstPath, 'dir') && bOverwrite
        if bVerbose
            warning('We merge folders now, be careful, this could go wrong when e.g. NIfTIs exist both unzipped & zipped');
        end
    end

    if isunix || ismac
        SrcPath = xASL_adm_UnixPath(SrcPath);
        DstPath = xASL_adm_UnixPath(DstPath);
		if exist(SrcPath, 'dir') && exist(DstPath, 'dir')
			% If we copy a directory to a directory in Linux, then this has to be done carefully, we want to copy the contents of SrcPath inside the 
			% DstPath. This has to be done using /*. Otherwise, it would have copied the entire SrcPath and put it inside DstPath
			if bOverwrite
				system(['cp -r -f ' SrcPath '/* ' DstPath]); % -n is short for --noclober
			else
				system(['cp -r ' SrcPath '/* ' DstPath]);
			end
        elseif exist(SrcPath, 'file') && exist(DstPath, 'dir')
            % If Srcpath is a file but DstPath is a folder, this goes
            % wrong, we would delete a folder when bOverwrite, so we throw
            % a warning and return
            warning('Trying to copy a file but the destination already exists as folder');
            fprintf('%s\n', 'Either delete the folder or append a file name, if the file should be copied into the folder');
        else
			% If one of SrcPath and DstPath is a file, then we copy normally
			if bOverwrite
				system(['cp -r -f ' SrcPath ' ' DstPath]); % -n is short for --noclober
			else
				system(['cp -r ' SrcPath ' ' DstPath]);
			end
		end
    else
        if exist(SrcPath, 'file') || exist(SrcPath, 'dir') % Check if source file exists
            if bOverwrite
                copyfile(SrcPath, DstPath, 'f');
            else
                copyfile(SrcPath, DstPath);
            end
        else
            warning('Source file does not exist...');
        end
    end
end


