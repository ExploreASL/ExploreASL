% Prepare EPAD QC QuickTest

%% THIS SCRIPT DELETES STUFF!!!!!!!

piet=pwd;
cd ../..
ExploreASL_Master('',0);
cd(piet);
Root = 'C:\Backup\ASL\EPAD\raw';

DirList = xASL_adm_GetFileList(Root, '.*', 'FPListRec', [0 Inf], true);

for iD=1:length(DirList)
    xASL_TrackProgress(iD, length(DirList));
    DcmList = xASL_adm_GetFileList(DirList{iD}, '.*\.dcm$', 'FPList', [0 Inf]);
    if ~isempty(DcmList)
        for iL=11:length(DcmList)
            xASL_delete(DcmList{iL});
        end
    end
end