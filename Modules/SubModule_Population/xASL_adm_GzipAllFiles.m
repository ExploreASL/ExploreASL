function xASL_adm_GzipAllFiles(ROOT, EXTERNAL, bFolder, bUseLinux)
%xASL_adm_GzipAllFiles Zip files or folders
%
% FORMAT: xASL_adm_GzipAllFiles(ROOT, bFolder, bUseLinux)
%
% INPUT:
%   ROOT        - the root study directory (REQUIRED)
%   EXTERNAL    - the path to the ExploreASL/External (OPTIONAL)
%   bFolder     - boolean, true for zipping complete folders rather than
%                 individual files (Unix only, OPTIONAL, DEFAULT=false)
%   bUseLinux   - boolean, true for using Unix's own filesystem and Gzip functionality.
%                 Has a much faster i/o than Matlab's i/o. Assumes that
%                 gzip is installed on Linux/MacOS (which is a fair
%                 assumption). (OPTIONAL, DEFAULT = true for Unix (Linux/Mac) and false
%                 otherwise (e.g. pc/Windows)
% OUTPUT: n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function zips NIfTI files or folders recursively and deletes
% the original file/folder after zipping.
%
% EXAMPLE: xASL_adm_GzipAllFiles('/MyStudy');
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


    %% ----------------------------------------------------
    %% 0) Admin
    if nargin<2 || isempty(EXTERNAL)
        EXTERNAL = '';
    end
    if nargin<3 || isempty(bFolder)
        bFolder = false;
    end
    if (nargin<4 || isempty(bUseLinux)) && (isunix || ismac)
        bUseLinux = true;
    elseif nargin<4 || isempty(bUseLinux)
        bUseLinux = false;
    end

    %% ----------------------------------------------------
    %% 1) Faster unix version, using OS file system
    if bUseLinux
        PathToSearch = xASL_adm_UnixPath(ROOT);
        system(['for i in `find ' PathToSearch '/* | grep -E \.nii$`; do gzip -1 -f -q -v "$i"; done'], '-echo');
        return;
    else
        %% 2) Otherwise use the multithreaded SuperGzip for Windows
        PathList = xASL_adm_GetFileList(ROOT, '^.*\.(nii)$', 'FPListRec', [0 Inf], false);
        fprintf('\n%s\n',['G-zZzZipping ' num2str(length(PathList)) ' files']);
        n_THREADS = min([str2num(getenv('NUMBER_OF_PROCESSORS')) length(PathList)]);
        if ~isempty(EXTERNAL)
            PathToSuperGzip = fullfile(EXTERNAL, 'SuperGZip', 'SuperGZip_Windows.exe');
            command = [PathToSuperGzip ' -p 0 -n ' num2str(n_THREADS) ' -v 1 ' ROOT ' *.nii'];
            exit_code = system(command);
        else
            exit_code = 1;
        end
        if exit_code == 0
            fprintf('Gzipping Niftis Successful...\n')
        else
            warning('An error occurred while using SuperGzip')
        end
    end

    %% ----------------------------------------------------
    %% 2) Slower Matlab version
    if bFolder
        PathList = xASL_adm_GetFileList(ROOT, '^.*$', 'FPList', [0 Inf], true);
        fprintf('\n%s\n',['G-zZzZipping ' num2str(length(PathList)) ' folders']);
%     else
%         PathList = xASL_adm_GetFileList(ROOT, '^.*\.(nii|wav|bmp)$', 'FPListRec', [0 Inf], false);
%         fprintf('\n%s\n',['G-zZzZipping ' num2str(length(PathList)) ' files']);
        
    end
    
%     fprintf('\n%s\n',['G-zZzZipping ' num2str(length(PathList)) ' files']);
%     fprintf('Progress:   ');
% 
%     for iList=1:length(PathList)
%         xASL_TrackProgress(iList,length(PathList));
%         
%         try
%             if bFolder
%                zip([PathList{iList} '.zip'], PathList{iList});
%                xASL_adm_DeleteFileList(PathList{iList}, '^.*$', true, [0 Inf]);
%                xASL_delete(PathList{iList});
%             else           
% 
%                 Fname = PathList{iList};
%                 [~, ~, Fext] = xASL_fileparts(Fname);
%                 if strcmp(Fext,'.nii') || strcmp(Fext,'.wav') || strcmp(Fext,'.bmp')
%                     if exist([Fname '.gz'],'file')
%                         delete([Fname '.gz']); 
%                     end
%                     gzip(Fname);
%                     delete(Fname);
%                 end
%             end
%         catch ME
%             warning('Something went wrong, continuing with next one');
%             fprintf('%s\n',ME.message);
%         end
%     end
%     fprintf('\n');

    
end