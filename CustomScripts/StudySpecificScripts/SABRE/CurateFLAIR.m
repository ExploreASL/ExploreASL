Root = '/Users/henk/Downloads/XNAT_FLAIR/cmp.slms.ucl.ac.uk/xnat/SABREv3';
FileList = xASL_adm_GetFileList(Root, '^.*\.dcm$', 'FPListRec' ,[0 Inf]);

for iFile=1:length(FileList)
    xASL_TrackProgress(iFile, length(FileList));
    ConvertDicomFolderStructure_CarefulSlow((fileparts(FileList{iFile})),0,1);
end

Root = '/Users/henk/ExploreASL/ASL/SABRE/raw';
FileList = xASL_adm_GetFileList(Root, '^.*\.dcm$', 'FPListRec' ,[0 Inf]);

for iFile=1:length(FileList)
    xASL_TrackProgress(iFile, length(FileList));
    if isempty(regexp(FileList{iFile},'FLAIR'))
        ConvertDicomFolderStructure_CarefulSlow((fileparts(FileList{iFile})),0,1);
    end
end