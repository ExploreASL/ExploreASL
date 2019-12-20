%% Redo Novice FLAIRs
% All FLAIRs & WMH exist in Longitudinal folder, with correct name
% So:

% 1) Overwrite FLAIR & WMH_SEGM, & delete lock files
dDir = 'C:\BackupWork\ASL\Novice\LongitudinalAnalysis';
oDir = 'C:\BackupWork\ASL\Novice\NoviceFLAIRbaseline';

SubjectList = xASL_adm_GetFileList(oDir,'^NOV\d{3}_FLAIR_biascorrected\.nii$','List');

for iSubject=1:length(SubjectList)
    SubjectID = SubjectList{iSubject}(1:6);
    
    NewFLAIR = fullfile(dDir,[SubjectID '_1'],'FLAIR.nii');
    NewWMH = fullfile(dDir,[SubjectID '_1'],'WMH_SEGM.nii');
    OldFLAIR = xASL_adm_GetFileList(oDir, ['^' SubjectID '_FLAIR_biascorrected\.nii$']);
    OldWMH = xASL_adm_GetFileList(oDir, ['^' SubjectID '_FLAIR_combined.*\.nii$']);
    
    LockDir = fullfile(dDir, 'lock', 'xASL_module_Structural', [SubjectID '_1'],'xASL_module_Structural');
    xASL_adm_DeleteFileList(LockDir, '.*', false, [0 Inf]);
    
    if length(OldFLAIR)~=1 || length(OldWMH)~=1
        error('CheckThis');
    else
        xASL_Copy(OldFLAIR{1}, NewFLAIR, true);
        xASL_Copy(OldWMH{1}, NewWMH, true);
    end
end


% 2) rerun the structural script & population script, as in the custom
% scripts folder

% 3) Recreate Figure 2

% 4) Zip all
% 5) Move all FLAIRs (baseline & followup) to Ldisk
% 6) move segmentations to follow-up FLAIR/Segmentations. "Provenance:
% these were manually segmented based on their registration of native space
% segmentation to follow up, by Malon vd ... and checked by Charles Majoie

%% 5) Copy FLAIRs & segmentations for FU
oDir = 'C:\BackupWork\ASL\Novice\LongitudinalAnalysis';
dDir = 'C:\BackupWork\ASL\Novice\FollowUpFLAIRs';

SubjectList = xASL_adm_GetFileList(oDir,'^NOV\d{3}_2$','List',[0 Inf],1);

for iSubject=1:length(SubjectList)
    FLAIRpathOrig = fullfile(oDir, SubjectList{iSubject}, 'FLAIR.nii');
    WMHpathOrig = fullfile(oDir, SubjectList{iSubject}, 'WMH_SEGM.nii');
    FLAIRpathDest = fullfile(dDir, ['FLAIR_' SubjectList{iSubject}(1:end-2) '.nii']);
    WMHpathDest = fullfile(dDir, ['WMH_SEGM_' SubjectList{iSubject}(1:end-2) '.nii']); 
    
    if xASL_exist(FLAIRpathOrig,'file')
        xASL_Copy(FLAIRpathOrig, FLAIRpathDest);
    end
    if xASL_exist(WMHpathOrig,'file')
        xASL_Copy(WMHpathOrig, WMHpathDest);
    end    
end
    
