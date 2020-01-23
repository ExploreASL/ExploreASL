% Concatenate NIfTIs

ASLdir = '/Users/henk/Downloads/sub01/ASL/ASL_2';
pathNewNifti = fullfile(ASLdir,'ASL4D.nii');
list = xASL_adm_GetFileList(ASLdir, '.*\.nii');

IM(:,:,:,1) = xASL_io_Nifti2Im(list{1});
IM(:,:,:,2) = xASL_io_Nifti2Im(list{5});
IM(:,:,:,3) = xASL_io_Nifti2Im(list{2});
IM(:,:,:,4) = xASL_io_Nifti2Im(list{6});
IM(:,:,:,5) = xASL_io_Nifti2Im(list{3});
IM(:,:,:,6) = xASL_io_Nifti2Im(list{7});
IM(:,:,:,7) = xASL_io_Nifti2Im(list{4});
IM(:,:,:,8) = xASL_io_Nifti2Im(list{8});

xASL_io_SaveNifti(list{1}, pathNewNifti, IM, [], 0);

PWI = IM(:,:,:,[1:2:end-1]) - IM(:,:,:,[2:2:end-0]);
dip_image(PWI)