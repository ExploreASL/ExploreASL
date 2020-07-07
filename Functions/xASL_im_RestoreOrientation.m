function xASL_im_RestoreOrientation(PathNIfTI)
%xASL_im_RestoreOrientation This function reverts the NIfTI header orientation matrix
% to the original orientation from the scanner/dcm2nii conversion
%
% FORMAT:       xASL_im_RestoreOrientation(PathNIfTI)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function reverts the NIfTI header orientation matrix
%               to the original orientation from the scanner/dcm2nii conversion.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL

if xASL_exist(PathNIfTI, 'file')
    nii = xASL_io_ReadNifti(PathNIfTI);
    nii.mat = nii.mat0;
    create(nii);
end
[Fpath, Ffile] = xASL_fileparts(PathNIfTI);
PathMat = fullfile(Fpath, [Ffile '.mat']);
xASL_delete(PathMat); % remove the orientation of other volumes when 4D as well

end
