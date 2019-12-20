function xASL_delete(InputPath, bFolderContent)
% Deletes the file given in the InputPath. For NIFTII file, delete also the GZ version
%
% FORMAT: xASL_delete(InputPath)
%
% INPUT:
%   InputPath           - path to the file to be deleted (REQUIRED)
%   bFolderContent      - TRUE for deleting all contents of a folder (if
%                         InputPath to folder) (OPTIONAL, DEFAULT = false)
%   
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Delete the file in the given path. If a NIFTI file with extension '.nii' or '.nii.gz' is given,
%              Then delete both the .nii and .nii.gz files.
% EXAMPLE: xASL_delete('/path/file.nii');
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright 2015-2019 ExploreASL
%
% 2019-09-29 HJ allow deleting folder including its contents

if nargin<2 || isempty(bFolderContent)
    bFolderContent = false;
elseif ischar(bFolderContent) || ~(bFolderContent==0 || bFolderContent==1)
    error(['Invalid parameter bFolderContent:' bFolderContent]);
end

if  xASL_exist(InputPath,'file')
	
	[Fpath, Ffile, Fext] = xASL_fileparts(InputPath);
	
	if ~isempty(strfind(Fext,'.nii'))
		% if we find a NIfTI, delete both .nii & .nii.gz
		InputNII = fullfile(Fpath,[Ffile '.nii']);
		InputGZ = fullfile(Fpath,[Ffile '.nii.gz']);
		
		if exist(InputNII,'file')
			delete(InputNII);
		end
		if exist(InputGZ ,'file')
			delete(InputGZ );
		end
	elseif exist(InputPath,'dir')
        % if we remove the full folder content:
        if bFolderContent
            xASL_adm_DeleteFileList(InputPath, '.*$', true, [0 Inf]); % remove all files in subdirectories
            DirList = xASL_adm_GetFileList(InputPath, '.*$', 'FPListRec', [0 Inf], true); % get list of folders
            for iDir=1:length(DirList) % delete all subfolders (assuming they are empty)
                xASL_delete(DirList{iDir});
            end
        end
            
        % first check whether the folder is empty % i.e. contains no files, empty subdirs we don't bother
        if bFolderContent || isempty(xASL_adm_GetFileList(InputPath, '.*$', 'FPListRec', [0 Inf], false))
            try 
                rmdir(InputPath, 's');
            catch
                warning([InputPath ' was empty but could not be deleted']);
            end
        else
            warning([InputPath ' is a non-empty directory']);
        end
    else
		delete(InputPath);
	end
	
end


end
