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
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

if xASL_exist(PathNIfTI, 'file')
    nii = xASL_io_ReadNifti(PathNIfTI);
    nii.mat = nii.mat0;
    create(nii);
end
[Fpath, Ffile] = xASL_fileparts(PathNIfTI);
PathMat = fullfile(Fpath, [Ffile '.mat']);
xASL_delete(PathMat); % remove the orientation of other volumes when 4D as well
end
