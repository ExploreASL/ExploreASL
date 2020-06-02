Rdir = '/Users/henk/ExploreASL/ASL/SABRE/Analysis2/Population';
FileList = xASL_adm_GetFileList(Rdir, '.*\.nii$', 'FPListRec', [0 Inf]);
for iFile=1:length(FileList)
    xASL_TrackProgress(iFile,length(FileList));
    xASL_adm_ZipFileNameHandling(FileList{iFile}, [], 1);
end