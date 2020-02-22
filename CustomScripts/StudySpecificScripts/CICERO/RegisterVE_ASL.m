%% Register VS_ASL

% 1) load DataPar Cicero

for iSubject=1:x.nSubjects
    SubjectDir = fullfile(x.D.ROOT, x.SUBJECTS{iSubject});
    VEdirs = xASL_adm_GetFileList(SubjectDir,'^VE_ASL_\d$','FPList',[0 Inf], true);
    
    PathT1 = fullfile(SubjectDir,'T1.nii');
    if ~xASL_exist(PathT1,'file')
        warning(['Not existing:' PathT1]);
        fprintf('\n\n\n');
        continue;
    end
    
    for iVE=1:length(VEdirs)
        PathMeanControl = fullfile(VEdirs{iVE},'mean_control.nii');
        if ~xASL_exist(PathMeanControl,'file')
            warning(['Not existing:' PathMeanControl]);
            fprintf('\n\n\n');
            continue;
        end
        xASL_adm_DeleteFileList(VEdirs{iVE}, '.*\.roi$', 1, [0 Inf]);
        % Find all NIfTIs except for the mean control
        OtherList = xASL_adm_GetFileList(VEdirs{iVE}, '.*(?<!mean_control)\.nii', 'FPListRec', [0 Inf]);
        
        % Register
        x.Quality = 1;
        xASL_spm_coreg(PathT1, PathMeanControl, OtherList, x);
    end
end