function xASL_Copy(SrcPath, DstPath, bOverwrite, bVerbose)
%xASL_Copy Copies a file or a directory, managing different OS, 
%.nii|.nii.gz & faster than Matlab implementation
%
% FORMAT: xASL_Copy(SrxASL_SysCopyath, DstPath[, bOverwrite, bVerbose)])
%
% INPUT:
%   SrcPath    - Source file or folder
%   DstPath    - Destination file or folder
%   bOverwrite - When true, overwrite if destination exists
%                When false, do nothing in that case (DEFAULT)
%   bVerbose   - If true then verbose output (DEFAULT=false);
%   
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Copies a file to a file or a directory to a directory. For a file, it zips it in the end if the destination path contains nii.gz. 
%              It also makes sure that only one of .nii and .nii.gz exists
%              in the destination directory. It is faster than the default
%              Matlab function, and runs on multiple OSes.
% 
%              NB: This function calls xASL_SysMove for the actual copying.
%              Run xASL_SysMove instead of xASL_Move if you don't want the
%              .nii|.nii.gz management/checking
%
% EXAMPLE: xASL_Copy('c:\User\path\file.nii', 'c:\User\path2\file.nii');             No overwriting
%          xASL_Copy('c:\User\path\file.nii.gz', 'c:\User\path2\file.nii.gz',true);  Overwriting and simple copying
%          xASL_Copy('c:\User\path\file.nii', 'c:\User\path2\file.nii.gz',true);     Overwriting and zipping
%          xASL_Copy('c:\User\path', 'c:\User\path2');                               Copying directories
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2019 ExploreASL

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

	% Obtains the paths and filenames of src and dest
    [~, SrcFile, SrcExt] = xASL_fileparts(SrcPath);
    [~, ~, DstExt] = xASL_fileparts(DstPath);

    %% Manage .nii|.nii.gz
    
	% If the destination is .nii.gz and source .nii, then zip the file in the end
    if ~exist(SrcPath,'dir') % does it only for file
        ZipInTheEnd = strcmp(DstExt, '.nii.gz') && strcmp(SrcExt, '.nii');

        % .gz compatibility ExploreASL
        [SrcPath, DstPath] = xASL_adm_ZipFileNameHandling(SrcPath, DstPath);
	else % not for whole directories
        ZipInTheEnd = false;
	end
	
	% If src is a file and destination a dir, then add the filename to the dir
	if ~exist(SrcPath,'dir') && exist(DstPath,'dir')
		DstPath = fullfile(DstPath,[SrcFile, SrcExt]);
    end
	
    %% Start copying
    
    % Checks if the destination file exists (as .nii or .nii.gz)
    if xASL_exist(DstPath,'file') 
        if bOverwrite && exist(DstPath,'dir')
			% Destination is a folder
            xASL_SysCopy(SrcPath, DstPath, bOverwrite, bVerbose); 
        elseif bOverwrite % destination is a file
            if bVerbose; fprintf('xASL_Copy: Overwriting %s\n', DstPath); end
            if exist(fileparts(DstPath),'dir') || isempty(fileparts(DstPath))
                xASL_delete(DstPath); % this will also delete .gz
            end
            xASL_adm_CreateDir(fileparts(DstPath)); % Create folder if doesnt exist
            xASL_SysCopy(SrcPath, DstPath, bOverwrite, bVerbose); 
        else
            if bVerbose; fprintf(2,'xASL_Copy: Skipping, destination exists: %s\n', DstPath); end
        end
	else
		% Destination folder or path does not exist - simply copy with given overwrite preferences
        xASL_adm_CreateDir(fileparts(DstPath)); % Create folder if doesnt exist
        xASL_SysCopy(SrcPath, DstPath, bOverwrite, bVerbose); 
        
		% Zip the destination file
        if ZipInTheEnd
            xASL_adm_GzipNifti(DstPath);

            if exist(DstPath,'file') && exist([DstPath '.gz'],'file')
                delete(DstPath);
            end
        end        
    end

end



