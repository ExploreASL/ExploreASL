%% Erode WM mask to small QC region where WM noise should be homogeneous & non-physiological

OriPath     = 'C:\ExploreASL\Maps\Atlases\DeepWM.nii';
OriIM       = xASL_io_Nifti2Im(OriPath);
QCpath      = 'C:\ExploreASL\Maps\Atlases\CentralWM_QC.nii';

NewIM                       = xASL_im_DilateErodeFull(OriIM,'erode',xASL_im_DilateErodeSphere(2));
NewIM                       = xASL_im_DilateErodeFull(NewIM,'erode',xASL_im_DilateErodeSphere(1));
NewIM(:,:,[1:67 76:end])    = 0; % mask out
NewIM(:,[1:62 83:end],:)    = 0;

xASL_io_SaveNifti(OriPath, QCpath, NewIM, [], 0);

% dip_image([OriIM NewIM])
