%% Pseudonymisation

% % Create key
% Pseudo = randn(300,1);
% Pseudo(:,2) = [1:300];
% Pseudo = sortrows(Pseudo,1);
% Pseudo(:,1) = [1:300];
% 
% for iSubject=1:length(Pseudo)
%     Key{iSubject,1} = ['NOV' sprintf('%03d',iSubject)];
%     Key{iSubject,2} = ['BlindFLAIR' sprintf('%03d',Pseudo(iSubject,2))];
% end
% 
% KeyPath = 'C:\BackupWork\ASL\Novice\Copy\FLAIR_intekening_okt2019\KeyBlindFLAIR.mat';
KeyPath = '/home/hjmutsaerts/lood_storage/divi/Projects/novice/MRI/FLAIR_intekening_okt2019/KeyBlindFLAIR.mat';

% save(KeyPath,'Key');

% clear

%% Move FLAIRs from NOV to pseudonymized
load(KeyPath,'-mat');

oDir = 'C:\BackupWork\ASL\Novice\Pseudonymization\flair_manualsegm';
dDir = 'C:\BackupWork\ASL\Novice\Pseudonymization\Pseudomyzed';

xASL_adm_CreateDir(dDir);
SubjList = xASL_adm_GetFileList(oDir,'^NOV\d{3}_(1|2)$','List',[0 Inf], true);

for iSubject=1:length(SubjList)
    xASL_TrackProgress(iSubject,length(SubjList));
    Subject = SubjList{iSubject}(1:end-2);
    TimePoint = SubjList{iSubject}(end);
    
    KeyIndex = find(cellfun(@(x) strcmp(x, Subject), Key(:,1)));
    
    DirOld = fullfile(oDir,SubjList{iSubject});
    DirNew = fullfile(dDir,[Key{KeyIndex,2} '_' TimePoint]);
    xASL_adm_CreateDir(DirNew);
    
    % FLAIR
    PathOld = fullfile(DirOld, 'rFLAIR.nii.gz');
    PathNew = fullfile(DirNew, [Key{KeyIndex,2} '_rFLAIR.nii.gz']);
    
    if xASL_exist(PathOld,'file')
        xASL_Copy(PathOld, PathNew, true);
    end
    % Segmentation
    FileSegm = xASL_adm_GetFileList(DirOld,'^NOV\d{3}_FLAIR_combined.*\.nii$', 'List');
    if ~isempty(FileSegm)
        PathOld = fullfile(DirOld, FileSegm{1});
        PathNew = fullfile(DirNew, [Key{KeyIndex,2} '_combined_label.nii.gz']);
        xASL_Copy(PathOld, PathNew, true);
    end
end

%% Move FLAIR & segmentations back (while backing up automatic WMH_SEGM)
load(KeyPath,'-mat');

oDir{2} = 'C:\BackupWork\ASL\Novice\Copy\FLAIR_intekening_okt2019\Pseudomyzed'; % follow up manual segmentations
oDir{1} = 'C:\BackupWork\ASL\Novice\NoviceFLAIRbaseline'; % baseline manual segmentations
dDir = 'C:\BackupWork\ASL\Novice\analysis';
NDir = 'C:\BackupWork\ASL\Novice\LongitudinalAnalysis';
xASL_adm_CreateDir(NDir);

SubjectList = xASL_adm_GetFileList(oDir{2},'^BlindFLAIR\d{3}_1$','List',[0 Inf], true);

for iSubject=1:length(SubjectList)
    xASL_TrackProgress(iSubject, length(SubjectList));
    % reconciliate name
    PseudoName = SubjectList{iSubject};
    IndexIs = find(cellfun(@(y) strcmp(y, PseudoName(1:end-2)), Key(:,2)));
    RealName{1} = [Key{IndexIs, 1} '_1'];
    RealName{2} = [Key{IndexIs, 1} '_2'];
    for TimePoint=1:2
        % move original files to longitudinal separate dir
        OldDir = fullfile(dDir, RealName{TimePoint});
        NewDir = fullfile(NDir, RealName{TimePoint});
        
        PathFLAIR = fullfile(NewDir, 'FLAIR.nii');
        PathWMH = fullfile(NewDir, 'WMH_SEGM.nii');
        PathFLAIRBackup = fullfile(NewDir, 'FLAIR_LST.nii');
        PathWMHBackup = fullfile(NewDir, 'WMH_LST.nii');
        
        if exist(OldDir,'dir')
            xASL_Move(OldDir, NewDir);
            % backup LST files
            xASL_Move(PathFLAIR,PathFLAIRBackup);
            xASL_Move(PathWMH,PathWMHBackup);
        end

        % move manual segmentations
        HadNoFLAIRorWMH = '';
        HadTooManyFLAIRorWMH = '';
        if TimePoint==1
            DirFLAIR = fullfile(oDir{1}, RealName{1}(1:end-2));
            DirWMH = DirFLAIR;
            PathFLAIRmanual = xASL_adm_GetFileList(DirFLAIR, '^NOV\d{3}_FLAIR.*biascorrected\.nii$');
            PathWMHmanual = xASL_adm_GetFileList(DirWMH, '^NOV\d{3}_FLAIR.*combined.*\.nii$');
        else
            DirFLAIR = fullfile(oDir{2}, [PseudoName(1:end-2) '_2']);
            DirWMH = fullfile(oDir{2}, PseudoName);
            PathFLAIRmanual = xASL_adm_GetFileList(DirFLAIR, ['^' PseudoName(1:end-2) '.*rFLAIR\.nii$']);
            PathWMHmanual = xASL_adm_GetFileList(DirWMH, ['^' PseudoName(1:end-2) '.*FU_segm\.nii$']);
        end
            
        if isempty(PathFLAIRmanual) || isempty(PathWMHmanual)
            HadNoFLAIRorWMH{end+1} = RealName{1};
        elseif length(PathFLAIRmanual)>1 || length(PathWMHmanual)>1
            HadTooManyFLAIRorWMH{end+1} = RealName{1};
        else
            xASL_Move(PathFLAIRmanual{1}, PathFLAIR);
            xASL_Move(PathWMHmanual{1}, PathWMH);
            % delete old folder
            if exist(DirFLAIR,'dir')
                xASL_adm_DeleteFileList(DirFLAIR, '.*', true, [0 Inf]);
                try
                    rmdir(DirFLAIR);
                end
            end
            if exist(DirWMH,'dir')
                xASL_adm_DeleteFileList(DirWMH, '.*', true, [0 Inf]);
                try
                    rmdir(DirWMH);
                end
            end
        end
    end
end
    
    
%% ExploreASL todo -> rerunning structural module by commenting everything else
% 1) Rerun 020_LinearReg_FLAIR2T1w only

% 2) add to provenance    
    
    
