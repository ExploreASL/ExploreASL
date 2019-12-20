%% Repair wrong T1 spatial definition/orientation

List            = {'AMC_23_2' 'AMC_32_1' 'AMC_34_2'};
OriginalT1w     = fullfile(x.D.ROOT,x.SUBJECTS{1},'T1_ORI.nii');

for iL=1:length(List)
    DFile       = fullfile(x.D.ROOT,List{iL},'T1.nii');
    IM          = xASL_io_Nifti2Im(DFile);
    xASL_io_SaveNifti(OriginalT1w, DFile, IM);
end

%% Repair wrong FLAIR spatial definition/orientation

List            = {'AMC_32_1' 'AMC_34_2'};
OriginalT1w     = fullfile(x.D.ROOT,x.SUBJECTS{1},'FLAIR_ORI.nii');

for iL=1:length(List)
    DFile       = fullfile(x.D.ROOT,List{iL},'FLAIR.nii');
    IM          = xASL_io_Nifti2Im(DFile);
    xASL_io_SaveNifti(OriginalT1w, DFile, IM);
end


%% RestoreOrientation

RestoreOrientation('C:\Backup\ASL\CP_Tavi\analysis\AMC_23_2\ASL_1\ASL4D.nii');

%% Repair wrong ASL spatial definition/orientation

List            = {'AMC_39_2'};
OriginalT1w     = fullfile(x.D.ROOT,x.SUBJECTS{3},'ASL_1','ASL4D.nii');

for iL=1:length(List)
    DFile       = fullfile(x.D.ROOT,List{iL},'ASL_1','ASL4D.nii');
    RestoreOrientation(DFile);
    OrientMat   = fullfile(x.D.ROOT,List{iL},'ASL_1','ASL4D.mat');
    if  exist(OrientMat,'file')
        delete(OrientMat);
    end

    IM          = xASL_io_Nifti2Im(DFile);
    xASL_io_SaveNifti(OriginalT1w, DFile, IM);

    OrientMat   = fullfile(x.D.ROOT,List{iL},'ASL_1','ASL4D.mat');
    if  exist(OrientMat,'file')
        delete(OrientMat);
    end
    RestoreOrientation(DFile);
end

%% RestoreSliceDirection
clear IMnew
for iS=1:size(IM,3)
    IMnew(:,:,iS,:)     = IM(:,:,size(IM,3)-iS+1,:);
end

clear IM
IM(:,:,1,:)                 = IMnew(:,:,end,:);
IM(:,:,2:size(IMnew,3),:)   = IMnew(:,:,1:end-1,:);

dip_image(IM)
