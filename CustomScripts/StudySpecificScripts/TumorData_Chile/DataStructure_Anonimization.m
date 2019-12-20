%% Anonimization folders

% Initiatlize ExploreASL without loading data, and specify your ROOT directory here

ROOT            = 'C:\Backup\ASL\Chile\raw';
AnalysisDir     = 'C:\Backup\ASL\Chile\analysis';
Dlist           = xASL_adm_GetFsList( ROOT,'^.*$',1)';

for iD=1:length(Dlist)
    OldName     = Dlist{iD};
    if  strcmp(OldName(1:3),'Sub')
        % Skip this
    else
        NewName     = ['Sub-' sprintf('%04d',iD)];
        OldDir      = fullfile(ROOT,Dlist{iD});
        NewDir      = fullfile(ROOT,NewName);
        xASL_Rename(OldDir,NewName);
        ID_key{iD,1}    = NewName;
        ID_key{iD,2}    = OldName;          
    end
end

ID_key          = ID_key(:,1:2);
ID_keyFilename  = fullfile(ROOT, 'ID_key.mat');
if  ~exist(ID_keyFilename,'file')
    save(ID_keyFilename,'ID_key');
end

% PS & XX & DS (Mac) files were deleted, these are not required

% -> hdr & img should be converted to nii
% use the command varargout=cat_io_img2nii(img,c,verb)

Flist   = xASL_adm_GetFileList( ROOT, '^2.*\.img$','FPListRec',[0 Inf]);
for iL=1:length(Flist)
    varargout   =   cat_io_img2nii(Flist{iL});
end

% Rename T1 niftis
Flist   = xASL_adm_GetFileList( ROOT, '^2.*T1.*\.nii$','FPListRec',[0 Inf]);
for iL=1:length(Flist)
    if  ~strcmp(Flist{iL},'T1.nii')
        xASL_Rename(Flist{iL},'T1.nii');
    end
end

% Rename voi to nii
Flist   = xASL_adm_GetFileList( ROOT, '^LESION.*\.voi$','FPListRec',[0 Inf]);
for iL=1:length(Flist)
    [Fpath Ffile Fext]  = fileparts(Flist{iL});
    NewFile     = [Ffile '.nii.gz'];
    xASL_Rename(Flist{iL},NewFile);
end

% Rename Lesion names
Flist   = xASL_adm_GetFileList( ROOT, '^LESION.*\.(nii|voi|nii.gz)$','FPListRec',[0 Inf]);
for iL=1:length(Flist)
    NewFile     = 'Lesion_T1_1.nii.gz';
    xASL_Rename(Flist{iL},NewFile);
end

% Remove useless files
Flist   = xASL_adm_GetFileList( ROOT, '^._.*$','FPListRec',[0 Inf]);
for iL=1:length(Flist)
    delete(Flist{iL});
end

% Get names from DICOMs -> this part takes very long, but it does work
Flist   = xASL_adm_GetFileList( ROOT, '^IM','FPListRec',[0 Inf]);
for iL=1:length(Flist)
    [Fpath Ffile Fext]  = fileparts(Flist{iL});
    ConvertDicomFolderStructure( Fpath );
end


% If a subject has both pCASL_1500ms and SOURCE__pCASL_1500ms, you should
% take the source one
Flist   = xASL_adm_GetFileList( ROOT, 'pCASL','FPListRec',[0 Inf]);

%% Create missing directories
Dlist   = xASL_adm_GetFsList(ROOT, '^Sub-\d{4}$',1);
for iD=1:length(Dlist)
    NewDir  = fullfile(AnalysisDir, Dlist{iD});
    xASL_adm_CreateDir(NewDir);
end


%% Copy the T1 NIfTIs from raw to analysis directory
Flist   = xASL_adm_GetFileList( ROOT, '^T1\.nii$','FPListRec',[0 Inf]);
for iL=1:length(Flist)
        try
            xASL_TrackProgress(iL,length(Flist));
            [Fpath Ffile Fext]  = fileparts(Flist{iL}(1:end-7));
            NewFile             = fullfile(AnalysisDir, Ffile, 'T1.nii');
            if  ~exist(NewFile,'file')
                xASL_Copy(Flist{iL},NewFile,1);
            end
        end
end

%% Copy the Lesion NIfTIs from raw to analysis directory
Flist   = xASL_adm_GetFileList( ROOT, '^Lesion.*\.(nii|nii\.gz)$','FPListRec',[0 Inf]);
for iL=1:length(Flist)
    try
        xASL_TrackProgress(iL,length(Flist));
        [Fpath Ffile Fext]  = xASL_fileparts(Flist{iL});
        FileName            = [Ffile Fext];
        [Fpath Ffile Fext]  = xASL_fileparts(Fpath);
        NewFile             = fullfile(AnalysisDir, Ffile, FileName);
        if  ~exist(NewFile,'file')
            xASL_Copy(Flist{iL},NewFile,1);
        end
    end
end
