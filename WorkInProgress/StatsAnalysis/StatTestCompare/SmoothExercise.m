% Copyright 2015-2024 ExploreASL (Works In Progress code)
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

%% Smooth exercise

addpath(genpath('c:\ASL_pipeline_HJ'));

temp = xASL_io_ReadNifti('E:\Backup\ASL_E\KCL_stats\INOX\smoothExercise\DARTEL_c1T1_20002.nii');

FwHm        = 8/1.5;
FwHm2SD     = (2*(2*reallog(2))^0.5);
SD          = FwHm/FwHm2SD;

%dip_image(temp.dat(:,:,:))
%dip_image(dip_arrtempIM)

%tempIM = gaussf ( temp.dat(:,:,:) , [SD SD SD], 'fir');
xASL_im_ndnanfilter(temp.dat(:,:,:),'gauss',[SD SD SD]*2.335,0);

ORIGINAL_NII    = 'E:\Backup\ASL_E\KCL_stats\INOX\smoothExercise\DARTEL_c1T1_20002.nii';
FILEPATH_NEW    = 'E:\Backup\ASL_E\KCL_stats\INOX\smoothExercise\DARTEL_c1T1_20002NEW.nii';
NEW_IMAGE       = tempIM;

%xASL_io_SaveNifti( ORIGINAL_NII , FILEPATH_NEW, dip_array(tempIM));
xASL_io_SaveNifti( ORIGINAL_NII , FILEPATH_NEW, tempIM);
