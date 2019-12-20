ExploreASL_Master('',0);

ROOT = 'C:\BackupWork\ASL\Dent';
RawFolder = fullfile(ROOT,'raw');
AnalysisFolder = fullfile(ROOT,'analysis');
SubjectList = xASL_adm_GetFileList(RawFolder,'^MCI-\d{4}$','FPListRec',[0 Inf], true);

%% Create foldernames
for iSubject=1:length(SubjectList)
    ConvertDicomFolderStructure_CarefulSlow(SubjectList{iSubject});
end
%% Replace spaces by underscore in foldernames
FolderList = xASL_adm_GetFileList(RawFolder,'^.*.$','FPListRec',[0 Inf], true);
for iFolder = 1:length(FolderList) 
    [Fpath, Fname, Fext] = xASL_fileparts(FolderList{iFolder});
    Fname = [Fname Fext];
    CorrectedName = xASL_adm_CorrectName(Fname,1,'_-\/'); 
    PathNew = fullfile(Fpath, CorrectedName);
    PathOld = fullfile(Fpath, Fname);
    if ~isequal(PathNew,PathOld)
		xASL_Move(PathOld, PathNew);
    end
end

%% Remove the XX & PS dicom files
xASL_adm_DeleteFileList(RawFolder, '^(PS|XX)_.*\.dcm$', true, [0 Inf]);

%% Remove empty folders
FolderList = xASL_adm_GetFileList(RawFolder,'^.*.$','FPListRec',[0 Inf], true);
for iFolder = 1:length(FolderList)
    if isempty(xASL_adm_GetFileList(FolderList{iFolder},'.*','FPListRec',[0 Inf], false))
        rmdir(FolderList{iFolder});
    end
end

%% Run import
imPar = ExploreASL_ImportConfig(ROOT);
%% Explore ASL Import script
bCopySingleDicoms = 0;
bUseDCMTK = false;
bCheckPermissions = false;
ExploreASL_Import(imPar, bCopySingleDicoms, bUseDCMTK, bCheckPermissions);

%% Fix M0
FileList = xASL_adm_GetFileList(AnalysisFolder, '^ASL4D\.nii$','FPListRec',[0 Inf]);
for iASL=1:length(FileList)
    xASL_io_SplitASL_M0(FileList{iASL},[1 2]); % first two volumes are M0
end

