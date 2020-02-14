%% Recreate Figure2 NOVICE but with identical flowfields for 2nd TP as for 1st TP


%% 1) Backup *_2/y_T1.nii & clone *_1/y_T1.nii > _2/y_T1.nii
for iSubj=1:x.nTimePointSubjects(1)
    xASL_TrackProgress(iSubj, x.nTimePointSubjects(1));
    y1 = fullfile(x.ROOT,x.TimePointSubjects{1}{iSubj},'y_T1.nii');
    y2 = fullfile(x.ROOT,x.TimePointSubjects{2}{iSubj},'y_T1.nii');
    y2Backup = fullfile(x.ROOT,x.TimePointSubjects{2}{iSubj},'y_T1_Backup.nii');

    xASL_Move(y2, y2Backup, 1);
    xASL_Copy(y1, y2, 1);
    
%     VALIDATION
%     IM1 = xASL_io_Nifti2Im(y1);
%     IM2 = xASL_io_Nifti2Im(y2);
%     CheckIM = IM1==IM2;
%     if min(CheckIM(:))~=1
%         warning(num2str(iSubj));
%     end
end

%% 2) Do the same for the WMH_SEGM.nii files
% Back them up, & replace their mat with the mat from TP1
for iSubj=1:x.nTimePointSubjects(1)
    xASL_TrackProgress(iSubj, x.nTimePointSubjects(1));
    WMH1 = fullfile(x.ROOT,x.TimePointSubjects{1}{iSubj},'WMH_SEGM.nii');
    WMH2 = fullfile(x.ROOT,x.TimePointSubjects{2}{iSubj},'WMH_SEGM.nii');
    WMH2Backup = fullfile(x.ROOT,x.TimePointSubjects{2}{iSubj},'WMH_SEGM_Backup.nii');

    xASL_Copy(WMH2, WMH2Backup, 1); % backup
    % replace mat
    nii1 = xASL_io_ReadNifti(WMH1);
    nii2 = xASL_io_ReadNifti(WMH2);
    nii2.mat = nii1.mat;
    nii2.mat0 = nii1.mat0;
    create(nii2);
    clear nii1 nii2
end
    

%% 2) Remove lock files for TP2 resampling only
LockFilesAre = {'080_Resample2StandardSpace' '100_VisualQC_Structural' '999_ready'};
for iSubj=1:x.nTimePointSubjects(1)
    xASL_TrackProgress(iSubj, x.nTimePointSubjects(1));
    LockFolder = fullfile(x.ROOT, 'lock', 'xASL_module_Structural', x.TimePointSubjects{2}{iSubj}, 'xASL_module_Structural');
    for iLock=1:length(LockFilesAre)
        LockPath = fullfile(LockFolder, [LockFilesAre{iLock} '.status']);
        xASL_delete(LockPath);
    end
end

%% 3) Rerun resampling for TP2 only, with structural xASL_module_Structural_Novice, only steps 
%% 4) rerun xASL_module_Population_NoviceLongitudinal

%% 5) restore the *_2/y_T1.nii from *_2/y_T1_Backup.nii
for iSubj=1:x.nTimePointSubjects(1)
    xASL_TrackProgress(iSubj, x.nTimePointSubjects(1));

    y2 = fullfile(x.ROOT,x.TimePointSubjects{2}{iSubj},'y_T1.nii');
    y2Backup = fullfile(x.ROOT,x.TimePointSubjects{2}{iSubj},'y_T1_Backup.nii');

    xASL_Move(y2Backup, y2, 1);
    
%     if xASL_exist(y2Backup,'file')
%         warning([num2str(iSubj) ' incorrectly copied']);
%     end
end

%% 6) restore the WMH
for iSubj=1:x.nTimePointSubjects(1)
    xASL_TrackProgress(iSubj, x.nTimePointSubjects(1));

    WMH2 = fullfile(x.ROOT,x.TimePointSubjects{2}{iSubj},'WMH_SEGM.nii');
    WMH2Backup = fullfile(x.ROOT,x.TimePointSubjects{2}{iSubj},'WMH_SEGM_Backup.nii');

    xASL_Move(WMH2Backup, WMH2, 1);
    
%     if xASL_exist(WMH2Backup,'file')
%         warning([num2str(iSubj) ' incorrectly copied']);
%     end    
end