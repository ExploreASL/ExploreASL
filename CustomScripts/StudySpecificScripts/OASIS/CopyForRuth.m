%% Copy CBF/PVC for Ruth Oliver

Odir = 'C:\BackupWork\ASL\OASIS';
Ddir = 'C:\BackupWork\ASL\OASIS\ForRuth';
AnalysisDirs = {'analysis_Siemens_Biograph_mMR_51010_syngo_MR_B18P' 'analysis_Siemens_Biograph_mMR_51010_syngo_MR_B20P' 'analysis_Siemens_TrioTim_35177_syngo_MR_B17' 'analysis_Siemens_TrioTim_35248_syngo_MR_B17'};

for iAnalysis=4 % 1:length(AnalysisDirs)
    OriDir = fullfile(Odir, AnalysisDirs{iAnalysis});
    DstDir = fullfile(Ddir, AnalysisDirs{iAnalysis});
    xASL_adm_CreateDir(DstDir);
    
    bMove = false; % Copy
    bDir = false; % search for files only
    bRecursive = true;
    bOverwrite = false;
    bVerbose = true;
    
    xASL_adm_CopyMoveFileList(OriDir, DstDir, '^(PVgm|PVwm|PVcsf|CBF|LabelingTerritories|MaskVascular)\.nii$', bMove, bDir, bRecursive, bOverwrite, bVerbose);
end

xASL_adm_GzipAllFiles('C:\BackupWork\ASL\OASIS\ForRuth');