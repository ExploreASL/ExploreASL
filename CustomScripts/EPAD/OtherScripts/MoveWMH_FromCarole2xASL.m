%% Move WMH_SEGM EPAD to analysis folder

ExploreASL_Master('',0);

iDir    = 'C:\Backup\ASL\EPAD\WMHcarole\WMH_Carole';
oDir    = 'C:\Backup\ASL\EPAD\analysis';

%% 1. Replace 'Correct_WS3WTC3WC1Lesion' with 'WMH_SEGM'
Flist   = xASL_adm_GetFileList(iDir,'^CorrectedCheckSCS_Correct_WS3WT3WC1Lesion.*\.(nii|nii\.gz)$');

for iL=1:length(Flist)
    iPath   = Flist{iL};
    [Fpath, Ffile, Fext]    = xASL_fileparts(Flist{iL});
    oPath                   = fullfile(Fpath, ['WMH_SEGM' Ffile(42:51) Fext]);
    if  exist(oPath,'file')
        oPath               = fullfile(Fpath, ['WMH_SEGM' Ffile(42:51) '_2' Fext]);
    end
    xASL_Move(iPath,oPath);
end

%% Move files into EPAD file structure
Dlist   = xASL_adm_GetFsList(oDir,'^\d{3}-\d{5}$',1);
Flist   = xASL_adm_GetFileList(iDir,'^WMH_SEGM_.*\.(nii|nii\.gz)$');

for iJ=1:length(Flist)
    [~, Ffile, ~]   = xASL_fileparts(Flist{iJ});
    SubjJ{iJ,1}   = Ffile(10:18);
end

NoSegmList = '';

for iL=2:length(Dlist)
    IndexN  = find(strcmp(SubjJ,Dlist{iL}));
    FileO   = fullfile(oDir,Dlist{iL},'WMH_SEGM.nii.gz');
    if ~isempty(IndexN)
        FileI   = Flist{IndexN};
        xASL_Move( FileI, FileO,1);
    end
    
    if ~exist(FileO,'file')
        NoSegmList{end+1,1}   = Dlist{iL};
    end
end
        
Flist   = xASL_adm_GetFileList(iDir,'^WMH_SEGM_.*\.(nii|nii\.gz)$');
for iL=1:length(Flist)
    [~, Ffile, ~]       = xASL_fileparts(Flist{iL});
    NoFLAIRlist{iL}     = Ffile(10:18);
end
NoFLAIRlist             = NoFLAIRlist';
