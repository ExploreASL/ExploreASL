DicomDirList = xASL_adm_GetFileList('/Users/henk/ExploreASL/ASL/AFIVASC_PedroVilela', '^DICOM$', 'FPListRec', [0 Inf], 1);
bUseDCMTK = false;
for iDir=1:length(DicomDirList)
    xASL_TrackProgress(iDir,length(DicomDirList));
    ConvertDicomFolderStructure_CarefulSlow(DicomDirList{iDir}, bUseDCMTK, 0);
end