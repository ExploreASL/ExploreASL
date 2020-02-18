%% Curate 7T Sikkel data
ExploreASL_Master('',0);

Rdir = '/home/hjmutsaerts/lood_storage/divi/Projects/sikkel7t/7T/';

Dlist = xASL_adm_GetFileList(Rdir,'^.*$', 'FPList', [0 Inf], true);

for iList=1:length(Dlist)
    xASL_TrackProgress(iList, length(Dlist));
    DICOMdir = fullfile(Dlist{iList}, 'DICOM');
    %% 1) Move metadata to separate dir
    MetaDir = fullfile(DICOMdir, 'MetadataDICOM');
    xASL_adm_CreateDir(MetaDir);
    MetaList = xASL_adm_GetFileList(DICOMdir, '^(PS|XX).*$','FPListRec',[0 Inf]);
    for iMeta=1:length(MetaList)
        [~, Ffile, Fext] = fileparts(MetaList{iMeta});
        xASL_Move(MetaList{iMeta}, fullfile(MetaDir,[Ffile Fext]),[],0);
    end
    
    %% 1.5) delete empty folders
    FolderList = xASL_adm_GetFileList(DICOMdir, '^.*$','FPListRec',[0 Inf], true);
    for iFolder=1:length(FolderList)
        if isempty(xASL_adm_GetFileList(FolderList{iFolder}, '^.*$','FPListRec',[0 Inf]))
            rmdir(FolderList{iFolder});
        end
    end
    
%     %% 2) add dicom suffixes
%     DICOMlist = xASL_adm_GetFileList(DICOMdir, '^(?<!PS)(?<!XX).*(?!\.dcm)$','FPListRec',[0 Inf]);
%     for iDICOM=1:length(DICOMlist)
%         [~, Ffile, Fext] = fileparts(DICOMlist{iDICOM});
%         xASL_Move(DICOMlist{iDICOM}, fullfile(DICOMdir,[Ffile '.dcm']));
%     end        
    %% 3) Sort in scantype folders
%     ConvertDicomFolderStructure_CarefulSlow(DICOMdir,1,0);
end

