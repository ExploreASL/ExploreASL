function [FilesList, FilesExeList, FoldersList] = xASL_adm_CheckPermissions(InputPath, FilesExecutable)
%xASL_adm_CheckPermissions Provide overview of permissions, try to change them if desired
%
% FORMAT: [FilesList, FilesExeList, FoldersList] = xASL_adm_CheckPermissions(InputPath[, FilesExecutable])
% 
% INPUT:
%   InputPath   - path to root folder containing files and folders to be checked for permissions (REQUIRED)
%   FilesExecutable - true if we want executable permissions for files (OPTIONAL, DEFAULT=false)
%
% OUTPUT:
%   FilesList     - cell structure table containing original file attributes
%   FilesExeList  - cell structure table containing original attributes from executable files
%   FoldersList   - cell structure table containing original folder attributes
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function does a recursive search through the root
%              folder & makes a list of the attributes of all files and folders.
%              It tries to reset the attributes to what we desire, which is by default:
%
%              - 664 for files (meaning only reading & writing for users & group, & read-only for others)
%              - 775 for folders (meaning reading, writing & opening for current user & current group, & for others only reading & opening)
%              
%              For executable files we also want 775.
%              Note that the permission to 'execute a folder' means opening them.
%              
%              - DataOK checks data permissions.
%              - ExeOK checks executable permissions.
%              - DataOK also includes executable permissions for folders.
%              - This runs recursively (but currently skips the contents of the root-folder) .
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_adm_CheckPermissions('/data/RAD/share/EPAD500_new/raw)
% __________________________________
% Copyright 2015-2020 ExploreASL

    %% Check permissions
    fprintf('%s\n', 'Checking permissions:   ');
    
    if nargin<2 || isempty(FilesExecutable)
        FilesExecutable = false; % by default, we don want executable permissions for files
    end
    
    if ~exist(InputPath, 'dir')
        error(['InputPath ' InputPath ' does not exist']);
	end
    
	% Read files and Dirs
    FileList = xASL_adm_GetFileList(InputPath, '.*', 'FPListRec', [0 Inf], false);
    DirList = xASL_adm_GetFileList(InputPath, '.*', 'FPListRec', [0 Inf], true);
    DirList{end+1} = InputPath;
	
	% Calculate number
    nFiles = length(FileList);
	nDirs  = length(DirList);

	% Reads the attributes for all files
    for iD=1:nFiles
        xASL_TrackProgress(iD, (nDirs + nFiles) * 2) % the *2 is to note that there may be some time needed to change permissions
        [~, attribN{iD}] = fileattrib(FileList{iD}); % try to read attributes
	end
    
	% Reads the attributes for all dirs
    for iD=1:nDirs
        xASL_TrackProgress(iD + nFiles, (nDirs + nFiles) * 2) % the *2 is to note that there may be some time needed to change permissions
        [~, attribN{iD+nFiles}] = fileattrib(DirList{iD});
    end
    
    % how do we call the permissions
    DataPermissions = {'UserRead' 'UserWrite' 'GroupRead' 'GroupWrite'};
    ExePermissions = {'UserExecute' 'GroupExecute'};

    FilesOK = true;
    FilesExeOK = true;
    FoldersOK = true;
    
    FilesList = '';
    FilesExeList = '';
    FoldersList = '';    
    
    % fill attrib
    if ~exist('attribN', 'var') || isempty(attribN)
        ExistAttribFile = false; % meaning that we couldnt read attributes above
    else
        ExistAttribFile = true;
    end
    
    if ~ExistAttribFile
        xASL_TrackProgress(2, 2);
        fprintf('\n');        
        warning(['Couldnt read attributes for ' InputPath ', skipping...']);
        return;
    else
        attrib = attribN{1};
        for iD=2:length(attribN)
            Ind1 = length(attrib)+1;
            Ind2 = length(attrib)+length(attribN{iD});
            attrib(Ind1:Ind2) = attribN{iD};
        end
    end
    

    % here we go through the attributes we found, & see if they comply with
    % what we wanted
    for iD=1:length(attrib)
        for iP=1:length(DataPermissions)
            if ~isnan(attrib(iD).(DataPermissions{iP})) % skip NaNs (e.g. on Windows filesystems)
                if ~attrib(iD).(DataPermissions{iP})
                    if attrib(iD).directory % for a folder
                        FoldersOK = false;
                        FoldersList{end+1,1} = attrib(iD).Name(length(InputPath)+1:end);
                    else
                        FilesOK = false;
                        FilesList{end+1,1} = attrib(iD).Name(length(InputPath)+1:end);
                    end
                end
            end
        end

        % here we do the same for the executable permissions
		for iP=1:length(ExePermissions)
			if ~isnan(attrib(iD).(ExePermissions{iP})) % skip NaNs (e.g. on Windows filesystems)
				if attrib(iD).directory % for a folder, we want to be able to execute (i.e. open) it
					if ~attrib(iD).(ExePermissions{iP})
						FoldersOK = false;
						FoldersList{end+1,1} = attrib(iD).Name(length(InputPath)+1:end);
					end
				else % if this is a file, then we check if we are expecting executable files
					if	attrib(iD).(ExePermissions{iP}) ~= FilesExecutable
						FilesExeOK = false;
						FilesExeList{end+1,1} = attrib(iD).Name(length(InputPath)+1:end);
					end
				end
			end
		end
    end
    
    %% Handle permissions for folders
    %  Here we try to change the folder permissions
    if ~FoldersOK 
        try
            xASL_TrackProgress(1, 2); % we are half way
            if isunix
                system(['find ' InputPath ' -type d -print0 | xargs -0 chmod 775']);
            elseif ispc
                system(['attrib -r ' InputPath ' /d /s']); % recursively removes read-only for directories on Windows
            end
            fprintf('%s\n', ['Folder permissions were reset recursively to 775 in ' InputPath]);
        catch
            fprintf('%s\n', 'For data we want to read and write, for users and groups (i.e. at least 660)');
            fprintf('%s\n', 'For executables/folders we want to read, write and execute for users and groups (i.e. at least 770)');
            error(['We do not have correct recursive permissions for folders in ' InputPath ', nor the ability to reset them']);
        end
    end
    
    %% Handle permissions for files
    %  Here we try to change the file permissions
    FileReset = false;
    if ~FilesExecutable % slight difference between wanting to be able to execute files as well
        ChmodPars = '664';
        if ~FilesOK
            FileReset = true;
        end
    else
        ChmodPars = '775';
        if ~FilesOK || ~FilesExeOK
            FileReset = true;
        end
    end
    
    if FileReset
        try
            xASL_TrackProgress(1, 2); % we are half way
            if isunix
                system(['find ' InputPath ' -type f -print0 | xargs -0 chmod ' ChmodPars]);
            elseif ispc
                system(['attrib -r ' InputPath '\*.* /s']); % recursively removes read-only for files on Windows
            end
            fprintf('%s\n', ['Files permissions were reset recursively to ' num2str(ChmodPars) 'in ' InputPath]);
        catch
            fprintf('%s\n', 'For data we want to read and write, for users and groups (i.e. at least 660)');
            fprintf('%s\n', 'For executables/folders we want to read, write and execute for users and groups (i.e. at least 770)');            
            error(['We do not have correct recursive permissions for files in ' InputPath ', nor the ability to reset them']);
        end
    end
    
    xASL_TrackProgress(2, 2);    
    fprintf('\n');
end

