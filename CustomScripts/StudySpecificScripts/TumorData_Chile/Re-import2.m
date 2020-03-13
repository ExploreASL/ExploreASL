%% Chile re-import

AnalysisDir     = 'C:\Backup\ASL\Chile\analysis';
RawDir          = 'C:\Backup\ASL\Chile\raw';
RegExp          = '^(Sub-\d{4}|Con-\d{3})$';

% 1) Move T1 & ASL NIfTIs to analysis -> from raw, then overwrite from ASL_Chile_patients

Dlist           = xASL_adm_GetFsList(RawDir,RegExp,1);

for iD=1:length(Dlist)
    xASL_TrackProgress(iD,length(Dlist));
    T1_File_ORI     = fullfile(RawDir,Dlist{iD},'T1.nii');
    T1_File_NEW     = fullfile(AnalysisDir,Dlist{iD},'T1.nii');
    if  xASL_exist(T1_File_ORI,'file')
        xASL_adm_CreateDir(fullfile(AnalysisDir,Dlist{iD}));
        xASL_Move(T1_File_ORI,T1_File_NEW,1);
    end

    segm_ORI        = fullfile(RawDir,Dlist{iD},'Lesion_T1_1.nii.gz');
    segm_NEW        = fullfile(AnalysisDir,Dlist{iD},'Lesion_T1_1.nii.gz');
    if  xASL_exist(segm_ORI,'file')
        xASL_adm_CreateDir(fullfile(AnalysisDir,Dlist{iD}));
        xASL_Move(segm_ORI,segm_NEW,1);
    end
end
