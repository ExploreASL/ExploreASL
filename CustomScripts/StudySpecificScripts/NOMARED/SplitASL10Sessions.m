%% ASL data splitting into 10 "sessions"
AnalysisDir = '/Volumes/MacUpSmall/NOMARED/analysis';
SubjectList = xASL_adm_GetFileList(AnalysisDir, '^NMRD\d{3}$', 'FPList', [0 Inf], 1);

for iSubject=1:length(SubjectList)
    xASL_TrackProgress(iSubject,length(SubjectList));
    for iPLD = 1:10
        % Manage folder
        DestDir = fullfile(SubjectList{iSubject}, ['ASL_' num2str(iPLD)]);
        DestPath = fullfile(DestDir, 'ASL4D.nii');
        xASL_adm_CreateDir(DestDir);
        
        % Load original files
        ControlPath = fullfile(SubjectList{iSubject}, 'ASL_1', ['ASL4D_' num2str(iPLD) '.nii']);
        LabelPath = fullfile(SubjectList{iSubject}, 'ASL_2', ['ASL4D_' num2str(iPLD) '.nii']);
        
        if xASL_exist(ControlPath,'file') && xASL_exist(LabelPath,'file')
            clear ASLim
            ASLim(:,:,:,1) = xASL_io_Nifti2Im(ControlPath);
            ASLim(:,:,:,2) = xASL_io_Nifti2Im(LabelPath);

            % Save new NIfTI
            xASL_io_SaveNifti(ControlPath, DestPath, ASLim, [], 0);
            % Delete old NIfTI
            xASL_delete(ControlPath);
            xASL_delete(LabelPath);
            % Move JSON
            PathJSONOrig = [ControlPath(1:end-4) '.json'];
            PathJSONDest = [DestPath(1:end-4) '.json'];
            xASL_Move(PathJSONOrig, PathJSONDest);
            xASL_delete([LabelPath(1:end-4) '.json']);
        end

    end
    
        % Delete parms.mat
        xASL_delete(fullfile(SubjectList{iSubject},'ASL_1','ASL4D_parms.mat'));
        xASL_delete(fullfile(SubjectList{iSubject},'ASL_2','ASL4D_parms.mat'));
        xASL_delete(fullfile(SubjectList{iSubject},'FLAIR_parms.mat'));
        xASL_delete(fullfile(SubjectList{iSubject},'T1_parms.mat'));
end