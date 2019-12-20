function xASL_im_RestoreOrientation( InputNII )
%xASL_im_RestoreOrientation This function reverts the NIfTI header orientation matrix
% to the original orientation from the scanner/dcm2nii conversion

if  xASL_exist(InputNII,'file')
    NII       = xASL_io_ReadNifti(InputNII);
    NII.mat   = NII.mat0;
    create(NII);
end

end
