function [IsEqualResolution] = xASL_im_CompareNIfTIResolutionXYZ(PathNIfTI1, PathNIfTI2)
%xASL_im_CompareNIfTIResolutionXYZ Compare 3D resolution between 2 NIfTIs
%
% FORMAT: [IsEqualResolution] = xASL_im_CompareNIfTIResolutionXYZ(PathNIfTI1, PathNIfTI2)
%
% INPUT:
%   PathNIfTI1  - Path to first NIfTI file (CHAR ARRAY, REQUIRED)
%   PathNIfTI2  - Path to second NIfTI file (CHAR ARRAY, REQUIRED)
%
% OUTPUT:
%   IsEqualResolution  - true if the 3D resolution is identical
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function checks whether the X, Y and Z resolution of a
%              NIfTI with any number of dimensions is equal. It rounds for 2 floating
%              points, for both NIfTIs, to ensure that the same precision is compared.
%
% EXAMPLE: IsEqualResolution = xASL_im_CompareNIfTIResolutionXYZ('/ASL/MyStudyName/Sub-001/c1T1.nii', '/ASL/MyStudyName/Sub-001/c1T1.nii');
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


nii1 = xASL_io_ReadNifti(PathNIfTI1);
nii2 = xASL_io_ReadNifti(PathNIfTI2);

res1 = nii1.hdr.pixdim(2:4);
res2 = nii2.hdr.pixdim(2:4);

IsEqualResolution = isequal(xASL_round(res1,2), xASL_round(res2,2));

end

