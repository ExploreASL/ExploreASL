function xASL_adm_GzipAllFiles(ROOT, bFolder, bUseLinux, pathExternal)
%xASL_adm_GzipAllFiles Zip files or folders
%
% FORMAT: xASL_adm_GzipAllFiles(ROOT, bFolder, bUseLinux, pathExternal)
%
% INPUT:
%   ROOT         - Root path of the folder structure you want to gzip recursively (REQUIRED)
%   bFolder      - Boolean, true for zipping complete folders rather than
%                  individual files (Unix only, OPTIONAL, DEFAULT=false)
%   bUseLinux    - Boolean, true for using Unix's own filesystem and Gzip functionality.
%                  Has a much faster i/o than Matlab's i/o. Assumes that
%                  gzip is installed on Linux/MacOS (which is a fair
%                  assumption). (OPTIONAL, DEFAULT = true for Unix (Linux/Mac) and false
%                  otherwise (e.g. pc/Windows)
%   pathExternal - Path to external Gzip tools like SuperGzip, used for Windows (.../ExploreASL/External) 
%                  (OPTIONAL, DEFAULT = [])
% OUTPUT: n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function zips NIfTI files or folders recursively and deletes
% the original file/folder after zipping.
%
% EXAMPLE: xASL_adm_GzipAllFiles(x);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


    %% ----------------------------------------------------
    %% 0) Admin
    if nargin<1 || isempty(ROOT)
        error('Please insert an root path...');
    end
    if nargin<2 || isempty(bFolder)
        bFolder = false;
    end
    if (nargin<3 || isempty(bUseLinux)) && (isunix || ismac)
        bUseLinux = true;
    elseif nargin<3 || isempty(bUseLinux)
        bUseLinux = false;
    end
    if nargin<4 || isempty(pathExternal)
        pathExternal = [];
    end

    exit_code = NaN;
    
    %% ----------------------------------------------------
    %% 1) Faster unix version, using OS file system
    if bUseLinux
        PathToSearch = xASL_adm_UnixPath(ROOT);
		if exist(PathToSearch,'dir') % Check if directory exists, otherwise a random current directory would be zipped
			oldPath = pwd;
			[exit_code, system_result] = system(['cd ' PathToSearch '; for i in `find * | grep -E \.nii$`; do gzip -1 -f -q -v "$i"; done']);
        end
    else
        %% 2) Otherwise use the multithreaded SuperGzip for Windows
        if ~isempty(pathExternal)
            PathList = xASL_adm_GetFileList(ROOT, '^.*\.(nii)$', 'FPListRec', [0 Inf], false);
            if ~isempty(PathList)
                fprintf('\n%s\n',['G-zZzZipping ' num2str(length(PathList)) ' files']);
                numCores = feature('numcores');
                % Get SuperGzip path
                PathToSuperGzip = fullfile(pathExternal, 'SuperGZip', 'SuperGZip_Windows.exe');
                % Define SuperGzip command
                command = [PathToSuperGzip ' -p 0 -n ' num2str(numCores) ' -v 1 ' ROOT ' *.nii'];
                % Run script
                [exit_code, system_result] = system(command);
                % Check if SuperGzip was successful

            end
        end
    end
    
    if isnan(exit_code)
        return; % nothing happened, nothing to zip
    elseif exit_code == 0
        fprintf('Gzipping of NIfTIs successful...\n');
        return;
    else
        warning('An error occurred trying to zip using the system: %s', system_result);
        fprintf('Trying slower Matlab version...\n');
    end    

    %% ----------------------------------------------------
    %% 2) Slower Matlab version
    if ~usejava('desktop')
        warning('No JavaVirtualMachine loaded, skipping zipping');
        return;
    end
    
    if bFolder
        PathList = xASL_adm_GetFileList(ROOT, '^.*$', 'FPList', [0 Inf], true);
        if ~isempty(PathList)
            fprintf('%s\n',['G-zZzZipping ' num2str(length(PathList)) ' folders']);
        end
    else
        PathList = xASL_adm_GetFileList(ROOT, '^.*\.(nii|wav|bmp)$', 'FPListRec', [0 Inf], false);
        if ~isempty(PathList)
            fprintf('%s\n',['G-zZzZipping ' num2str(length(PathList)) ' files']);
        end
    end
    
    % Iterate over the individual files/folders
    if ~isempty(PathList)
        fprintf('Progress:   ');
        
        for iList=1:length(PathList)
            xASL_TrackProgress(iList,length(PathList));

            try
                if bFolder
                   zip([PathList{iList} '.zip'], PathList{iList});
                   xASL_adm_DeleteFileList(PathList{iList}, '^.*$', true, [0 Inf]);
                   xASL_delete(PathList{iList});
                else           

                    Fname = PathList{iList};
                    [~, ~, Fext] = xASL_fileparts(Fname);
                    if strcmp(Fext,'.nii') || strcmp(Fext,'.wav') || strcmp(Fext,'.bmp')
                        if exist([Fname '.gz'],'file')
                            delete([Fname '.gz']); 
                        end
                        gzip(Fname);
                        delete(Fname);
                    end
                end
            catch ME
                warning('Something went wrong, continuing with next one...');
                fprintf('%s\n',ME.message);
            end
        end
        fprintf('\n');
    end

    
end