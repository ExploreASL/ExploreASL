%% Zeroing out imaging data for BIDS example datasets

ExploreASL_Initialize;

Rdir = 'C:\Users\henkj\Google Drive\WorkGeneral\StudiesProjects\BIDS_ASL\TestDataSet_ImportBIDS\raw';
Flist = xASL_adm_GetFileList(Rdir, '.*\.(nii|nii\.gz)$', 'FPListRec', [0 Inf], false);

for iL=1:length(Flist)
    xASL_TrackProgress(iL,length(Flist));
    tN=xASL_io_ReadNifti(Flist{iL});
    NewIm=zeros(size(tN.dat),'uint8');
    xASL_io_SaveNifti(Flist{iL},Flist{iL},NewIm,[],true);
end
