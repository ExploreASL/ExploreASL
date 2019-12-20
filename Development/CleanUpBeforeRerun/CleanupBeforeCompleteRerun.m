function CleanupBeforeCompleteRerun(InputPath)
%CleanupBeforeCompleteRerun Reverts previous run, deleting derivatives, while keeping raw data intact

% InputPath = 'C:\Backup\ASL\Novice\analysis';

%% 1) Delete derivatives

Subdirs = {'Population' 'lock'};

for iS=1:length(Subdirs)
    Subdir = fullfile(InputPath, Subdirs{iS});
    if exist(Subdir, 'dir')
        rmdir(Subdir, 's');
    end
end

Ftype = {'ples.*\.(nii|nii\.gz)' 'c(1|2)T1\.(nii|nii\.gz)' 'j_T1\.(nii|nii\.gz)' 'rT1\.(nii|nii\.gz)' 'WMH_SEGM\.(nii|nii\.gz)' 'WMH_SEGM_CleanUp\.(nii|nii\.gz)' 'y_T1\.(nii|nii\.gz)' '.*\.(log|pdf|txt|json|csv|tsv|ps)'};
Ftype2 = {'T1_seg8\.mat' 'ASL_module.log' 'despiked_ASL4D.nii' 'ASL4D.mat' 'mean_PWI_beforeMoCo.nii' 'mean_control_beforeMoCo.nii' 'despiked_ASL4D.mat' 'FoV.nii' 'mean_control.nii' 'mmean_control.nii' 'rp_ASL4D_BeforeSpikeExclusion.txt' 'SD_control_beforeMoCo.nii' 'SD_PWI_beforeMoCo.nii' 'slice_gradient.nii' 'SNR_control_beforeMoCo.nii' 'SNR_PWI_beforeMoCo.nii' 'mean_PWI_Clipped.nii' 'PWI.nii' 'mean_PWI_Clipped_sn.mat' 'rgrey.nii' 'rASL4D.nii' 'mask_ICV.nii' 'TempFilter_ASL4D.nii' 'TempFilter_despiked_ASL4D.nii'};
IndN = length(Ftype);
Ind2N = length(Ftype2);
Ftype(IndN+1:IndN+Ind2N) = Ftype2;
clear Ftype2

for iL=1:length(Ftype)
    xASL_TrackProgress(iL, length(Ftype));
    Flist = xASL_adm_GetFileList(InputPath, Ftype{iL}, 'FPListRec', [0 Inf]);
    for iL2=1:length(Flist)
        xASL_delete(Flist{iL2});
    end
end

%% 2) Rename backupped/original files
OriTypeList = {'FLAIR' 'T1'};
for iO=1:length(OriTypeList)
    OriList = xASL_adm_GetFileList(InputPath, ['^' OriTypeList{iO} '_ORI\.(nii|nii\.gz)$'], 'FPListRec', [0 Inf]);
    for iL=1:length(OriList)
        xASL_TrackProgress(iL, length(OriList));
        OriPath = OriList{iL};
        [Fpath, Ffile, Fext] = xASL_fileparts(OriPath);
        NonPath = fullfile(Fpath, [OriTypeList{iO} '.nii']);
        if xASL_exist(NonPath, 'file') && xASL_exist(OriPath, 'file')
            xASL_delete(NonPath);
            xASL_Move(OriPath, NonPath);
        end
    end
end
    
fprintf('\n');

end
