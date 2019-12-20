function EPAD_Manage_PE_Niftis(AnalysisDir)
%EPAD_Manage_PE_Niftis Make sure all PE niftis contain only a single 3D volume

SubjList = xASL_adm_GetFsList(AnalysisDir,'^\d{3}EPAD\d{4}$', true, [], [], [0 Inf]);
fprintf('%s','Managing PE NIfTIs:  0%');

for iSubj=1:length(SubjList)
    xASL_TrackProgress(iSubj,length(SubjList));
    CurrDir = fullfile(AnalysisDir, SubjList{iSubj});
    Flist = xASL_adm_GetFileList(CurrDir, '.*_(Norm|Rev)PE.*\.(nii|nii\.gz)$', 'FPListRec', [0 Inf], false);
    for iL=1:length(Flist)
        tN = xASL_io_Nifti2Im(Flist{iL});
        tN = tN(:,:,:,1);
        xASL_io_SaveNifti(Flist{iL}, Flist{iL}, tN, [], 0);
    end
end

end
