ExploreASL_Master('',0);
%% 
Rdir = '/Users/henk/Downloads/SABRE/cmp.slms.ucl.ac.uk/xnat/SABREv3';
SubjList = xASL_adm_GetFileList(Rdir, ['.*\d{5}.*'], 'FPList', [0 Inf], true);

for iSubj=1:length(SubjList)
    xASL_TrackProgress(iSubj, length(SubjList));
    List1 = xASL_adm_GetFileList(SubjList{iSubj}, '.*', 'FPList', [0 Inf], true);
    if length(List1)~=1
        warning(num2str(iSubj));
        fprintf('\n\n\n');
        continue;
    end
    ScanList = xASL_adm_GetFileList(List1{1}, '.*', 'FPList', [0 Inf], true);
    for iScan=1:length(ScanList)
        DICOMdir = fullfile(ScanList{iScan},'DICOM');
        ConvertDicomFolderStructure_CarefulSlow(DICOMdir,0, 0);
    end
end

%% Remove last two characters
Rdir = '/Users/henk/Downloads/SABRE/cmp.slms.ucl.ac.uk/xnat/SABREv3';
SubjList = xASL_adm_GetFileList(Rdir, ['\d{5}\d*(_|)(P|E|I)'], 'FPList', [0 Inf], true);

for iSubj=1:length(SubjList)
    xASL_TrackProgress(iSubj,length(SubjList));
    [~, Ffile] = fileparts(SubjList{iSubj});
    StartInd = regexp(Ffile,'(_|)(P|E|I)');
    if ~isempty(StartInd)
        NewFile{iSubj,1} = Ffile;
        NewFile{iSubj,2} = Ffile(1:StartInd-1);
        Character{iSubj,1} = NewFile{iSubj,2};
        Character{iSubj,2} = Ffile(end);
        
        PathNew = fullfile(Rdir, NewFile{iSubj,2});
        xASL_Move(SubjList{iSubj}, PathNew);
    end
end

SavePath = '/Users/henk/Downloads/SABRE/analysis/Character.mat';
save(SavePath,'Character');

%% Fix the M0 scan
AnalysisDir = '/s4ever/radG/RAD/share/SABRE/Analysis2';
fprintf('Generating folder list\n');
Dlist = xASL_adm_GetFileList(AnalysisDir,'\d*','FPList', [0 Inf], true);

fprintf('Fixing the M0 scans:   ')

for iDir=1:length(Dlist)
    xASL_TrackProgress(iDir,length(Dlist));
    PathM0 = fullfile(Dlist{iDir},'ASL_1','M0.nii');
    PathM0Backup = fullfile(Dlist{iDir},'ASL_1','M0_Backup.nii');
    
    if xASL_exist(PathM0Backup,'file')
        xASL_Copy(PathM0Backup, PathM0, true);
    elseif ~xASL_exist(PathM0,'file')
        continue;
    else
        xASL_Copy(PathM0, PathM0Backup);
    end
    
    IM = xASL_io_Nifti2Im(PathM0);
    if length(size(IM))~=4 || max(size(IM)~=[80 80 20 3])
        warning(['Size M0 differed: ' PathM0]);
    else
        try
            FullIM = reshape(IM,[size(IM,1) size(IM,2) size(IM,3)*size(IM,4)]);
            clear IM;
            nVol=3;
            for ii=1:nVol
                IM(:,:,:,ii) = FullIM(:,:,[ii:nVol:end-(nVol-ii)]);
            end
            IM = xASL_stat_MeanNan(IM,4);
            xASL_io_SaveNifti(PathM0, PathM0, IM, [], 0);
        catch ME
            warning(ME.message);
        end
    end
end
fprintf('\n');

%% Fix the ASL scan
% AnalysisDir = '/s4ever/radG/RAD/share/SABRE/analysis';
% Dlist = xASL_adm_GetFileList(AnalysisDir,'\d*','FPList', [0 Inf], true);

fprintf('Fixing the ASL scans:   ')

for iDir=1:length(Dlist)
    xASL_TrackProgress(iDir,length(Dlist));
    PathASL = fullfile(Dlist{iDir},'ASL_1','ASL4D.nii');
    PathASLBackup = fullfile(Dlist{iDir},'ASL_1','ASL4D_Backup.nii');
    
    if xASL_exist(PathASLBackup,'file')
        xASL_Copy(PathASLBackup, PathASL, true);
    elseif ~xASL_exist(PathASL,'file')
        continue;
    else
        xASL_Copy(PathASL, PathASLBackup);
    end
    
    IM = xASL_io_Nifti2Im(PathASL);
    if length(size(IM))<4
        warning(['Dim issue:' PathASL]);
        fprintf('\n\n\n');
    else
        if max(size(IM)~=[80 80 40 35])
            warning(['Size ASL differed: ' PathASL]);
            fprintf('\n\n\n');
        end
        FullIM = reshape(IM,[size(IM,1) size(IM,2) size(IM,3)*size(IM,4)]);
        nVol = size(IM,3)*size(IM,4)/20;

        clear IM;
        for ii=1:nVol
            IM(:,:,:,ii) = FullIM(:,:,[ii:nVol:end-(nVol-ii)]);
        end
        FullIM = reshape(IM,[size(IM,1) size(IM,2) size(IM,3)/2 size(IM,4)*2]);
        IM = zeros([size(IM,1) size(IM,2) size(IM,3) size(IM,4)]);
        IM(:,:,[1:2:end-1],:) = FullIM(:,:,:,1:end/2);
        IM(:,:,[2:2:end-0],:) = FullIM(:,:,:,end/2+1:end);

        xASL_io_SaveNifti(PathASL, PathASL, IM, [], 0);
    end
end
fprintf('\n');

%% Move all subjects without ASL to exclusion folder
AnalysisDir = '/Users/henk/ExploreASL/ASL/SABRE/analysis';
ExclusionDir = '/Users/henk/ExploreASL/ASL/SABRE/ExcludedNoASL';
xASL_adm_CreateDir(ExclusionDir);
Dlist = xASL_adm_GetFileList(AnalysisDir,'\d*','FPList', [0 Inf], true);

for iDir=1:length(Dlist)
    xASL_TrackProgress(iDir,length(Dlist));
    PathASL = fullfile(Dlist{iDir},'ASL_1','ASL4D.nii');
    [~, Ffile] = xASL_fileparts(Dlist{iDir});
    List{iDir,1} = Ffile;
    if ~xASL_exist(PathASL,'file')
        xASL_Move(Dlist{iDir}, fullfile(ExclusionDir, Ffile));
        List{iDir,2} = 0;
    else
        List{iDir,2} = 1;
    end
end
        
        
