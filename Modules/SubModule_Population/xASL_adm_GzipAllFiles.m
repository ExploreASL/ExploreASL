function xASL_adm_GzipAllFiles(ROOT, bFolder, bUseLinux)
%xASL_adm_GzipAllFiles This function zips all .nii
% files recursively into nii.gz.
% Note that the linux option cannot deal with spaces in filenames or
% directories

    if nargin<2 || isempty(bFolder)
        bFolder = false;
    end
    if nargin<3 || isempty(bUseLinux) && (isunix || ismac)
        bUseLinux = true;
    elseif nargin<3 || isempty(bUseLinux)
        bUseLinux = false;
    end

    if bUseLinux
        PathToSearch = xASL_adm_UnixPath(ROOT);
        system(['for i in `find ' PathToSearch '/* | grep -E \.nii$`; do gzip -1 -f -q -v "$i"; done'], '-echo');
        return;
    end
            
    
    if bFolder
        PathList = xASL_adm_GetFileList(ROOT, '^.*$', 'FPList', [0 Inf], true);
        fprintf('\n%s\n',['G-zZzZipping ' num2str(length(PathList)) ' folders']);
    else
        PathList = xASL_adm_GetFileList(ROOT, '^.*\.(nii|wav|bmp)$', 'FPListRec', [0 Inf], false);
        fprintf('\n%s\n',['G-zZzZipping ' num2str(length(PathList)) ' files']);
    end
    
    fprintf('\n%s\n',['G-zZzZipping ' num2str(length(PathList)) ' files']);
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
            warning('Something went wrong, continuing with next one');
            fprintf('%s\n',ME.message);
        end
    end
    fprintf('\n');

end


