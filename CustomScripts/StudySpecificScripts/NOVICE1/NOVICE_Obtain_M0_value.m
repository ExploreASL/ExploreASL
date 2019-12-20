%% Get M0-value NOVICE

loadDIR     = 'E:\M0\M0\M0\1134_M0_eveline\20140110_1.3.46.670589.11.42151.5.0.8340.2014011018411468006\MR\01201_NOVICE_NEW';

loadNII     = fullfile( loadDIR, '00001.nii');
loadMAT     = fullfile( loadDIR, 'dcm_values.mat');
temp        = xASL_io_ReadNifti( loadNII );
load(loadMAT);
im          = temp.dat(:,:,:,:)./(dcm_values{4} * dcm_values{5});

%
% threshold   = 3* std(im(:));
% dip_image(xASL_im_rotate([im;100000*(single(im>threshold))],90))
%%%%%%%

clear
loadDIR     = 'E:\M0\M0\M0\01212_M0_veronica\20140110_1.3.46.670589.11.42151.5.0.8340.2014011018232082000\MR\01201_NOVICE_NEW';

loadNII     = fullfile( loadDIR, '00001.nii');
loadMAT     = fullfile( loadDIR, 'dcm_values.mat');
temp        = xASL_io_ReadNifti( loadNII );
load(loadMAT);
im          = temp.dat(:,:,:,:)./(dcm_values{4} * dcm_values{5});

dip_image(xASL_im_rotate(im,90))
%%%%%%%%
clear
loadDIR     = 'E:\M0\M0\M0\1213_M0_jordy\20140110_1.3.46.670589.11.42151.5.0.8340.2014011019330054024\MR\01201_NOVICE_NEW';

loadNII     = fullfile( loadDIR, '00001.nii');
loadMAT     = fullfile( loadDIR, 'dcm_values.mat');
temp        = xASL_io_ReadNifti( loadNII );
load(loadMAT);
im          = temp.dat(:,:,:,:)./(dcm_values{4} * dcm_values{5});

dip_image(xASL_im_rotate(im,90))

%%%%%%%
clear
loadDIR     = 'E:\M0\M0\M0\12345_Sophie_Cohen\20140108_1.3.46.670589.11.42151.5.0.2800.2014010817243334000\MR\00901_NOVICE_NEW_16CH';

loadNII     = fullfile( loadDIR, '00001.nii');
loadMAT     = fullfile( loadDIR, 'dcm_values.mat');
temp        = xASL_io_ReadNifti( loadNII );
load(loadMAT);
im          = temp.dat(:,:,:,:)./(dcm_values{4} * dcm_values{5});

dip_image(xASL_im_rotate(im,90))


%%%%%%%
clear
loadDIR     = 'E:\M0\M0\M0\123465_Carline_Tacke\20140108_1.3.46.670589.11.42151.5.0.2800.2014010817501190006\MR\01001_NOVICE_NEW_16CH';

loadNII     = fullfile( loadDIR, '00001.nii');
loadMAT     = fullfile( loadDIR, 'dcm_values.mat');
temp        = xASL_io_ReadNifti( loadNII );
load(loadMAT);
im          = temp.dat(:,:,:,:)./(dcm_values{4} * dcm_values{5});

dip_image(xASL_im_rotate(im,90))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIND/AgeIV/Epod
clear
loadDIR     = 'E:\M0\M0\M0\1134_M0_eveline\20140110_1.3.46.670589.11.42151.5.0.8340.2014011018411468006\MR\00901_FIND_NEW';

loadNII     = fullfile( loadDIR, '00001.nii');
loadMAT     = fullfile( loadDIR, 'dcm_values.mat');
temp        = xASL_io_ReadNifti( loadNII );
load(loadMAT);
im          = temp.dat(:,:,:,:)./(dcm_values{4} * dcm_values{5});

dip_image(xASL_im_rotate(im,90))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Epod
clear
loadDIR     = 'E:\M0\M0\M0\1134_M0_eveline\20140110_1.3.46.670589.11.42151.5.0.8340.2014011018411468006\MR\00801_EPOD_NEW';

loadNII     = fullfile( loadDIR, '00001.nii');
loadMAT     = fullfile( loadDIR, 'dcm_values.mat');
temp        = xASL_io_ReadNifti( loadNII );
load(loadMAT);
im          = temp.dat(:,:,:,:)./(dcm_values{4} * dcm_values{5});

dip_image(xASL_im_rotate(im,90))
