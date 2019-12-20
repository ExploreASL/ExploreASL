%% Fix NIfTIs
RootIn = 'C:\BackupWork\ASL\OASIS\OASIS3_xASL';
DirList = xASL_adm_GetFileList(RootIn, '^OAS\d*_\d*$', 'List', [0 Inf], true);

%% Remove any NifTI duplicates from dcm2nii
for iDir=1:length(DirList)
    xASL_TrackProgress(iDir, length(DirList));
    ASLList = xASL_adm_GetFileList(fullfile(RootIn, DirList{iDir}), '^ASL_\d*$', 'List', [0 Inf], true);
    for iASL = 1:length(ASLList)
        ASLdir = fullfile(RootIn, DirList{iDir},['ASL_' num2str(iASL)]);
        ASLnii = sort(xASL_adm_GetFileList(ASLdir, '^ASL4D.*\.nii$','FPList',[0 Inf]));
        ASLjson = sort(xASL_adm_GetFileList(ASLdir, '^ASL4D.*\.json$','FPList',[0 Inf]));
        
        if length(ASLnii)>1 && length(ASLjson)>1
            LengthList = max(length(ASLnii),length(ASLjson));
            for iList=2:LengthList
                xASL_delete(ASLjson{iList});
                xASL_delete(ASLnii{iList});
            end
        end
    end
end

%% Separate the M0
for iDir=1:length(DirList)
    xASL_TrackProgress(iDir, length(DirList));
    ASLList = xASL_adm_GetFileList(fullfile(RootIn, DirList{iDir}), '^ASL_\d*$', 'List', [0 Inf], true);
    for iASL = 1:length(ASLList)
        ASLdir = fullfile(RootIn, DirList{iDir},['ASL_' num2str(iASL)]);
        ASLnii = sort(xASL_adm_GetFileList(ASLdir, '^ASL4D\.nii$','FPList',[0 Inf]));
        xASL_io_SplitASL_M0(ASLnii{iASL},1);
    end
end