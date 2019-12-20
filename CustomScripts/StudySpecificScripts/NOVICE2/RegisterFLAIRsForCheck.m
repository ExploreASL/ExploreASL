%% Register FLAIRs Novice1 manual segmentation & Novice2 automatic segmentation

ManualDir = 'C:\BackupWork\ASL\Novice\Flair_WMH_check';
AnalysisDir = 'C:\BackupWork\ASL\Novice\analysis';

SubjList = xASL_adm_GetFileList(ManualDir,'^NOV\d{3}$','List',[0 Inf],true);

MissingFields = '';

for iSubj=1:length(SubjList)
    OriDir = fullfile(ManualDir, SubjList{iSubj});
    NewDir = fullfile(AnalysisDir, [SubjList{iSubj} '_1']);
    FLAIRori = xASL_adm_GetFileList(OriDir, [SubjList{iSubj} '_(FLAIR_|)biascorrect(ed|ion)\.nii']);
    WMHori = xASL_adm_GetFileList(OriDir, [SubjList{iSubj} '_(FLAIR_|)combined(-|)label\.nii']);
    
    if ~(length(FLAIRori)==1) || ~(length(WMHori)==1)
        MissingFields{end+1,1} = SubjList{iSubj};
    else
        
        % find AnalysisDir
        FLAIRnew = fullfile(AnalysisDir, [SubjList{iSubj} '_1'], 'FLAIR.nii');
        WMHnew = fullfile(AnalysisDir, [SubjList{iSubj} '_1'], 'WMH_SEGM.nii');
        
        if ~xASL_exist(FLAIRnew,'file') || ~xASL_exist(WMHnew,'file')
            MissingFields{end+1,1} = SubjList{iSubj};
        else
            % copy
            Move_WMHnew = fullfile(OriDir, 'WMH_SEGM.nii');
            if xASL_exist(Move_WMHnew,'file')
                % skip
            else
                % register
                xASL_spm_coreg(FLAIRnew, FLAIRori{1}, {WMHori{1}});
                % copy
                xASL_Copy(WMHnew, Move_WMHnew, true);
            end
        end
    end
end
        
MissingFields = MissingFields';
        
        
        % copy
        