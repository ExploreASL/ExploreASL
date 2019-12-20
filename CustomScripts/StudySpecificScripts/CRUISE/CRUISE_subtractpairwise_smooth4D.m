%% CRUISE subtract pairwise & smooth 4D ASL file

ASL_File= fullfile( x.D.ROOT, x.SUBJECTS{1}, x.SESSIONS{1}, 'rtemp_despiked_ASL4D.nii');
ASL_im  = xASL_io_ReadNifti(ASL_File);
ASL_im  = ASL_im.dat(:,:,:,:);
[ control_im label_im]  = Check_control_label_order( ASL_im );

% Paired subtraction
ASL_im                  = control_im - label_im;

% Create 4D smoothing kernel
FWHM_image      = 12; % mm smoothing kernel for the 3D image
FWHM_temp       = 15; % n subtractions smoothing kernel for temporal dimension

FwHm2SD         = (2*(2*reallog(2))^0.5);

FwHm_image      = 8/1.5; % to divide by voxel-size for number of voxels
SD_image        = FwHm_image/FwHm2SD;

FwHm_temp       = FWHM_temp;
SD_temp         = FwHm_temp/FwHm2SD;

ASL_im_smoothed = dip_array(gaussf ( ASL_im , [SD_image SD_image SD_image SD_temp], 'fir'));

% dip_image(ASL_im_smoothed(:,:,57,:))
clear ASL_im control_im label_im
FILEPATH_NEW    = fullfile( x.D.ROOT, x.SUBJECTS{1}, x.SESSIONS{1}, 'rtemp_despiked_ASL4D_smoothed4D.nii');
xASL_io_SaveNifti( ASL_File, FILEPATH_NEW, ASL_im_smoothed, size(ASL_im_smoothed,4), 16, size(ASL_im_smoothed,5) )
