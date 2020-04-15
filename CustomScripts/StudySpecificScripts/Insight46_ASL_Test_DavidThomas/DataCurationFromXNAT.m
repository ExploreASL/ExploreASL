x = ExploreASL_Master('',0);

ROOT = '/Users/henk/Downloads/Insight46';

Flist = xASL_adm_GetFileList(ROOT,'.*\.zip','FPListRec');
for iList=1:length(Flist)
    xASL_TrackProgress(iList,length(Flist));
    cd(fileparts(Flist{iList}));
    unzip(Flist{iList});
    xASL_delete(Flist{iList});
end

DirList = xASL_adm_GetFileList(RawDir,'.*','FPList',[0 Inf], 1);
for iDir=1:3
    ConvertDicomFolderStructure_CarefulSlow(DirList{iDir}, 0, 1);
end

DirList = xASL_adm_GetFileList(ROOT,'.*','List',[0 Inf], 1);
RawDir = fullfile(ROOT,'raw');
xASL_adm_CreateDir(RawDir);
for iDir=1:3
    xASL_Move(fullfile(ROOT,DirList{iDir}) , fullfile(ROOT,'raw',DirList{iDir}) );
end
    

ExploreASL_Import(ExploreASL_ImportConfig(ROOT));

% Put a DataPar file in, then run ExploreASL