Flist = xASL_adm_GetFileList('/home/hjmutsaerts/scratch/Data','^FLAIR.*run.*\.(json|nii|nii\.gz|mat)$', 'FPListRec');

for iList=1:length(Flist)
    xASL_TrackProgress(iList, length(Flist));
    PathOri = Flist{iList};
    [Fpath, Ffile, Fext] = fileparts(PathOri);
    if strcmp(Fext,'.mat')
        xASL_delete(PathOri);
    else
        PathDst = fullfile(Fpath, ['FLAIR' Fext]);
        if ~xASL_exist(PathDst,'file')
            xASL_Move(PathOri,PathDst);
        else
            xASL_delete(PathOri);
        end
    end
end
        
Flist = xASL_adm_GetFileList('/home/hjmutsaerts/scratch/Data','^T1.*run.*\.(json|nii|nii\.gz|mat)$', 'FPListRec');

for iList=1:length(Flist)
    xASL_TrackProgress(iList, length(Flist));
    PathOri = Flist{iList};
    [Fpath, Ffile, Fext] = fileparts(PathOri);
    if strcmp(Fext,'.mat')
        xASL_delete(PathOri);
    else
        PathDst = fullfile(Fpath, ['T1' Fext]);
        if ~xASL_exist(PathDst,'file')
            xASL_Move(PathOri,PathDst);
        else
            xASL_delete(PathOri);
        end
    end
end
