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
KeyPath = '/home/hjmutsaerts/lood_storage/divi/Projects/novice/MRI/proc/novice2/FLAIR/Pseudomyzed/KeyBlindFLAIR.mat';
% KeyPath = '/home/hjmutsaerts/lood_storage/divi/Projects/novice/MRI/FLAIR_intekening_okt2019/KeyBlindFLAIR.mat';

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


%% Deblind FLAIR & segmentations
load(KeyPath,'-mat');

Rdir = '/home/hjmutsaerts/lood_storage/divi/Projects/novice/MRI/proc/novice2/FLAIR/Pseudomyzed';
SubjectList = xASL_adm_GetFileList(Rdir,'^BlindFLAIR\d{3}_(1|2)$','List',[0 Inf], true);

for iSubject=1:length(SubjectList)
    xASL_TrackProgress(iSubject, length(SubjectList));
    % reconciliate name
    IndexIs = find(cellfun(@(y) strcmp(y, SubjectList{iSubject}(1:end-2)), Key(:,2)));
    RealName = [Key{IndexIs, 1} SubjectList{iSubject}(end-1:end)];
    
    DirOld = fullfile(Rdir, SubjectList{iSubject});
    DirNew = fullfile(Rdir, RealName);
    xASL_Move(DirOld, DirNew);
end

xASL_adm_DeleteFileList(Rdir, '.*\.(itksnap|png|txt)$', true, [0 Inf]);

%% Delete _2 dirs (we assume that all scans are in alignment with the first timepoint, these scans were only there for reference)
Rdir = '/home/hjmutsaerts/lood_storage/divi/Projects/novice/MRI/proc/novice2/FLAIR';
SubjectList = xASL_adm_GetFileList(Rdir,'^NOV\d{3}_2$','FPList',[0 Inf], true);

for iSubject=1:length(SubjectList)
    xASL_TrackProgress(iSubject, length(SubjectList));
    xASL_adm_DeleteFileList(SubjectList{iSubject}, '.*', true, [0 Inf]);
    system(['rmdir ' SubjectList{iSubject}]);
end

%% Rename dirs to remove _1 suffix
Rdir = '/home/hjmutsaerts/lood_storage/divi/Projects/novice/MRI/proc/novice2/FLAIR';
SubjectList = xASL_adm_GetFileList(Rdir,'^NOV\d{3}_1$','FPList',[0 Inf], true);

for iSubject=1:length(SubjectList)
    xASL_TrackProgress(iSubject, length(SubjectList));
    xASL_Move(SubjectList{iSubject}, SubjectList{iSubject}(1:end-2));
end

%% Remove existing manual WMH & FLAIR, and manage to have 1 FLAIR at least
Ddir = '/home/hjmutsaerts/lood_storage/divi/Projects/novice/MRI/ExploreASLoutput/LongitudinalAnalysis';
DirList = xASL_adm_GetFileList(Ddir,'^NOV\d{3}_(1|2)$','FPList',[0 Inf], true);

for iDir=1:length(DirList)
    xASL_TrackProgress(iDir, length(DirList));
    
    % remove WMH_SEGM.nii (this keeps WMH_SEGM_LST.nii)
    WMHpath = fullfile(DirList{iDir}, 'WMH_SEGM.nii');
    xASL_delete(WMHpath);

    % if exist FLAIR_LST.nii & FLAIR.nii, move FLAIR_LST to FLAIR.nii
    FLAIRpath = fullfile(DirList{iDir}, 'FLAIR.nii');
    FLAIR_LSTpath = fullfile(DirList{iDir}, 'FLAIR_LST.nii');
    if xASL_exist(FLAIR_LSTpath, 'file')
        xASL_Move(FLAIR_LSTpath, FLAIRpath, true);
    end
    
    % if no FLAIR, give warning
    if ~xASL_exist(FLAIRpath,'file')
        warning(FLAIRpath);
        fprintf('\n\n\n');
    end
    % if WMH, give warning
    if xASL_exist(WMHpath,'file')
        warning(WMHpath);
        fprintf('\n\n\n');
    end    
end

%% Do the same for the non-longitudinal data
Ddir = '/home/hjmutsaerts/lood_storage/divi/Projects/novice/MRI/ExploreASLoutput/NonLongitudinal';
DirList = xASL_adm_GetFileList(Ddir,'^NOV\d{3}_(1|2)$','FPList',[0 Inf], true);

for iDir=1:length(DirList)
    xASL_TrackProgress(iDir, length(DirList));
    
    WMHpath = fullfile(DirList{iDir}, 'WMH_SEGM.nii');
    WMH_LSTpath = fullfile(DirList{iDir}, 'WMH_SEGM_LST.nii');
    FLAIRpath = fullfile(DirList{iDir}, 'FLAIR.nii');
    FLAIR_LSTpath = fullfile(DirList{iDir}, 'FLAIR_LST.nii');    
    
    % manage WMH
    if xASL_exist(WMH_LSTpath,'file') && xASL_exist(WMHpath,'file')
        % if both manual & LST segmentations exist, remove the manual (old
        % wrong one)
        xASL_delete(WMHpath);
    elseif xASL_exist(WMHpath,'file')
        % if only the LST one exist, back this up as such
        xASL_Move(WMHpath, WMH_LSTpath);
    end
        
    % if FLAIR_LST exists, move to FLAIR.nii
    if xASL_exist(FLAIR_LSTpath, 'file')
        xASL_Move(FLAIR_LSTpath, FLAIRpath, true);
    end
    
    % if no FLAIR, give warning
    if ~xASL_exist(FLAIRpath,'file')
        warning(FLAIRpath);
        fprintf('\n\n\n');
    end
end


%% Register FLAIRs & move WMHs if exist. if not existing, create empty WMH_segm
Ddir = '/home/hjmutsaerts/lood_storage/divi/Projects/novice/MRI/ExploreASLoutput/LongitudinalAnalysis';
Odir = '/home/hjmutsaerts/lood_storage/divi/Projects/novice/MRI/proc/novice2/FLAIR';

DirList = xASL_adm_GetFileList(Odir,'^NOV\d{3}$','List',[0 Inf], true);

for iDir=1:length(DirList)
    xASL_TrackProgress(iDir, length(DirList));
    
    % Define paths
    DdirCurrent1 = fullfile(Ddir, [DirList{iDir} '_1']);
    DdirCurrent2 = fullfile(Ddir, [DirList{iDir} '_2']);
    Path_FLAIR_Dest1 = fullfile(DdirCurrent1, 'FLAIR.nii');
    Path_WMH_Dest1 = fullfile(DdirCurrent1, 'WMH_SEGM.nii');
    Path_FLAIR_Dest2 = fullfile(DdirCurrent2, 'FLAIR.nii');
    Path_WMH_Dest2 = fullfile(DdirCurrent2, 'WMH_SEGM.nii');    
    
    if ~xASL_exist(Path_WMH_Dest1,'file') || ~xASL_exist(Path_WMH_Dest2,'file')
        % skip this iteration if already done
        OdirCurrent = fullfile(Odir, DirList{iDir});
        Path_CombiSegm = xASL_adm_GetFileList(OdirCurrent, '.*FU_segm.*\.nii$', 'FPList', [0 Inf]);
        Path_WrongSegm = xASL_adm_GetFileList(OdirCurrent, '.*combined_label\.nii$', 'FPList', [0 Inf]);
        Path_FLAIR = xASL_adm_GetFileList(OdirCurrent, '.*rFLAIR.*\.nii$', 'FPList', [0 Inf]);

        if length(Path_CombiSegm)==1 && length(Path_WrongSegm)==1 
            xASL_delete(Path_WrongSegm{1});
        end
        
        if ~(length(Path_CombiSegm)==1 && length(Path_FLAIR)==1)
            warning(['Couldnt find FLAIR & WMH for n=' DirList{iDir}]);
            fprintf('\n\n\n');
        elseif length(xASL_adm_GetFileList(OdirCurrent, '^.*$'))~=2
            warning(['Length wasnt 2 for n=' DirList{iDir}]);
            fprintf('\n\n\n');
        else
            % Define temporary paths
            [Fpath, Ffile] = xASL_fileparts(Path_FLAIR);
            Path_FLAIRtemp = fullfile(Fpath, [Ffile '_temp.nii']);
            [Fpath, Ffile] = xASL_fileparts(Path_CombiSegm{1});
            Path_CombiSegmTemp = fullfile(Fpath, [Ffile '_temp.nii']);    
            % Delete temporary paths
            xASL_delete(Path_FLAIRtemp);
            xASL_delete(Path_CombiSegmTemp);              
            
            % register FLAIR for timepoint 1 (should it already be)
            xASL_spm_coreg(Path_FLAIR_Dest1, Path_FLAIR{1}, {Path_CombiSegm{1}});
            % if CombiSegm==1 || 3, this is a baseline lesion
            tIM = xASL_io_Nifti2Im(Path_CombiSegm{1});
            if sum(tIM(:)<0)~=0 || sum(tIM(:)>3)~=0
                warning(['incorrect labels present in ' Path_CombiSegm{1}]);
            else
                tIM = tIM==1 | tIM==3;
                xASL_io_SaveNifti(Path_CombiSegm{1}, Path_WMH_Dest1, tIM);
            end
            
            % do the same for timepoint 2, but with temporary copy
            xASL_Copy(Path_FLAIR{1}, Path_FLAIRtemp, true);
            xASL_Copy(Path_CombiSegm{1}, Path_CombiSegmTemp, true);
            
            % register FLAIR for timepoint 2
            xASL_spm_coreg(Path_FLAIR_Dest2, Path_FLAIRtemp, {Path_CombiSegmTemp});
            % if CombiSegm==1 || 3, this is a baseline lesion
            tIM = xASL_io_Nifti2Im(Path_CombiSegmTemp);
            if sum(tIM(:)<0)~=0 || sum(tIM(:)>3)~=0
                warning(['incorrect labels present in ' Path_CombiSegmTemp]);
            else
                tIM = tIM>0;
                xASL_io_SaveNifti(Path_CombiSegmTemp, Path_WMH_Dest2, tIM);
            end
            xASL_delete(Path_FLAIRtemp);
            xASL_delete(Path_CombiSegmTemp);
        end
    end
end

%% Check for missing WMHs in the longitudinal dataset
Ddir = '/home/hjmutsaerts/lood_storage/divi/Projects/novice/MRI/ExploreASLoutput/LongitudinalAnalysis';
DirList = xASL_adm_GetFileList(Ddir,'^NOV\d{3}_(1|2)$','List',[0 Inf], true);

for iDir=1:length(DirList)
    xASL_TrackProgress(iDir, length(DirList));
    
    Path_WMH = fullfile(Ddir, DirList{iDir}, 'WMH_SEGM.nii');
    
    if ~xASL_exist(Path_WMH,'file')
        warning([Path_WMH ' didnt exist']);
        fprintf('\n\n\n');
    end
end
        
        
%% ExploreASL -> rerunning resampling, get tissue volume & visual QC

ExploreASL_Master('/home/hjmutsaerts/lood_storage/divi/Projects/novice/MRI/ExploreASLoutput/LongitudinalAnalysis/DATA_PAR.json',1,1,1,1,1);
    
% then move all other files from non-longitudinal into longitudinal folder
