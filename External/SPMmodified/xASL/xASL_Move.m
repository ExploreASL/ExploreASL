function xASL_Move(SrcPath, DstPath, bOverwrite, bVerbose)
%xASL_Move Move a file or a directory. Manage nii/gz & different OSes
%
% FORMAT: xASL_Move(SrcPath, DstPath[, bOverwrite, bVerbose])
%
% INPUT:
%   SrcPath    - Source file or folder (REQUIRED)
%   DstPath    - Destination file or folder (REQUIRED)
%   bOverwrite - When true, overwrite if destination exists
%                When false, do nothing in that case (OPTIONAL, DEFAULT = FALSE)
%   bVerbose   - Verbose mode (OPTIONAL, DEFAULT = TRUE)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Moves a file to a file, a file to a directory, or a directory to a directory. It keeps the initial extensions, no unzipping or zipping
%              after the move. But it makes sure that only one of .nii and .nii.gz exists in the destination directory.
%              Bypass inefficient matlab stuff on linux and windows, but
%              can only move on the same file system.
% 
%              NB: This function calls xASL_SysMove for the actual moving.
%              Run xASL_SysMove instead of xASL_Move if you don't want the
%              .nii|.nii.gz management/checking
% 
% EXAMPLE: xASL_Move('c:\User\path\file.nii', 'c:\User\path2\file.nii');             No overwriting
%          xASL_Move('c:\User\path\file.nii.gz', 'c:\User\path2\file.nii.gz',true);  Overwriting and simple copying
%          xASL_Move('c:\User\path\file.nii', 'c:\User\path2\file.nii.gz',true);     Overwriting and zipping
%          xASL_Move('c:\User\path', 'c:\User\path2',[],false);                      Moving directories, no overwrite, verbose off
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright (C) 2015-2021 ExploreASL

    %% Manage input arguments
	if nargin<4 || isempty(bVerbose)
		bVerbose = true;
    end
    if nargin<3 || isempty(bOverwrite)
        bOverwrite = false;
    end    
    if nargin<2
        error('Please provide both SrcPath & DstPath');
    elseif strcmp(SrcPath,DstPath)
        if bVerbose
            warning('SrcPath & DstPath were equal, no xASL_Move action required');
            fprintf('%s\n', SrcPath);
        end
        return;
    end
	% Check if source exists
	if ~xASL_exist(SrcPath,'file')
		error(['Source doesnt exist: ' SrcPath]);
	end
	
	% Obtains the paths and filenames of src and dest
    [~, SrcFile, SrcExt] = xASL_fileparts(SrcPath);
    [~, ~, DstExt] = xASL_fileparts(DstPath);
	
    %% Manage .nii|.nii.gz
    
	% If the destination is .nii.gz and source .nii, then zip the file in the end
	if ~exist(SrcPath,'dir') % does it only for file
        bZipInTheEnd = strcmp(DstExt, '.nii.gz') && strcmp(SrcExt, '.nii');

        % .gz compatibility ExploreASL
        [SrcPath, DstPath] = xASL_adm_ZipFileNameHandling(SrcPath, DstPath);
	else % not for whole directories
        bZipInTheEnd = false;
	end
	
    % If src is a file and destination a dir, then add the filename to the dir
	if exist(DstPath,'dir') && ~exist(SrcPath,'dir')
		DstPath = fullfile(DstPath, [SrcFile, SrcExt]);
	end
	
	% Checks if the destination file exists (as .nii or .nii.gz)
	if ~exist(fileparts(DstPath),'dir')
        xASL_adm_CreateDir(fileparts(DstPath));
    end
    
    %% Start moving
    if xASL_exist(DstPath,'file') % this will also check .gz
		if bOverwrite
			if bVerbose; fprintf('xASL_Move: overwriting %s\n', DstPath); end

			% When destination is a file, then remove both .nii and .nii.gz
			if ~exist(DstPath,'dir') && xASL_exist(DstPath,'file')
				xASL_delete(DstPath);
            end
            xASL_adm_CreateDir(fileparts(DstPath)); % Create folder if doesnt exist
			xASL_SysMove(SrcPath, DstPath, bOverwrite);
		else
			if bVerbose; fprintf(2,'xASL_Move: skip: destination exists: %s\n', DstPath); end
			bZipInTheEnd = false;
		end
    else
        xASL_adm_CreateDir(fileparts(DstPath)); % Create folder if doesnt exist
		xASL_SysMove(SrcPath, DstPath, bOverwrite);
	end
    
	% Zip the destination file
	if bZipInTheEnd
		xASL_adm_GzipNifti(DstPath);
		
		if exist(DstPath,'file') && exist([DstPath '.gz'],'file')
			delete(DstPath);
		end
	end
	
end
