function xASL_adm_GzipAllFiles(ROOT)
%xASL_adm_GzipAllFiles This function zips all temporary
% files.
%
% By default, the low priority files are zipped,
% also, the high priority files can be zipped (to be made default later?)
%
% Low priority files are those that are never used for re-processing, only sometimes checked
% visually
%
% High priority files are those that can be used for re-processing

    
    PathList    = xASL_adm_GetFileList(ROOT, '^.*\.(nii|wav|bmp)$', 'FPListRec', [0 Inf]);
    fprintf('\n%s\n',['G-zZzZipping ' num2str(length(PathList)) ' files']);
    fprintf('Progress:   ');

    for iP=1:length(PathList)
        xASL_TrackProgress(iP,length(PathList));
        Fname = PathList{iP};
        [~, ~, Fext] = xASL_fileparts(Fname);
        if strcmp(Fext,'.nii')
            try
                if exist([Fname '.gz'],'file')
                    delete([Fname '.gz']); 
                end
                gzip(Fname);
                delete(Fname);
            catch ME
                warning('Something went wrong, continuing with next one');
                fprintf('%s\n',ME.message);
            end
        end
    end
    fprintf('\n');

end


