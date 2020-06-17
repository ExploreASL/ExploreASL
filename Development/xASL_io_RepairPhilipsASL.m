function xASL_io_RepairPhilipsASL(InPath)
%xASL_im_RepairPhilipsASL Performs correction of Philips slice-time order
%from dcm2nii. Make sure to backup the ASL file first

    tIM     = xASL_io_Nifti2Im(InPath);
    cIM     = zeros(size(tIM));

    sliceHalf = size(tIM,3)/2;
    timeHalf = size(tIM,4)/2;
    cIM(:,:,1:2:end,1:2:end) = tIM(:,:,1:sliceHalf,1:timeHalf);
    cIM(:,:,1:2:end,2:2:end) = tIM(:,:,(sliceHalf+1):end,1:timeHalf);

    cIM(:,:,2:2:end,1:2:end) = tIM(:,:,1:sliceHalf,(timeHalf+1):end);
    cIM(:,:,2:2:end,2:2:end) = tIM(:,:,(sliceHalf+1):end,(timeHalf+1):end);

    xASL_io_SaveNifti(InPath,InPath,cIM,[],0);

end
