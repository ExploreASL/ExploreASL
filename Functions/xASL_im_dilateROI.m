function xASL_im_dilateROI(PathIn, PathTemp)
%xASL_im_dilateROI ...
%
% FORMAT:       xASL_im_dilateROI(PathIn, PathTemp)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  ...
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL

if nargin<2
    PathTemp = PathIn;
end

IM      = xASL_io_Nifti2Im(PathIn);
IM      = xASL_im_ConvertMap2Mask(IM);
HDR     = xASL_io_ReadNifti(PathIn);
pixdim  = prod(HDR.hdr.pixdim(2:4));
IMvol   = sum(IM(:)).*pixdim;
while  IMvol<40
       IM = xASL_im_DilateErodeFull(IM,'dilate',xASL_im_DilateErodeSphere(2));
       IMvol   = sum(IM(:)).*pixdim;
end

xASL_io_SaveNifti(PathIn, PathTemp, IM, [], false);

% For all regions smaller than 40 mm^3, dilate it until it is at least 40 mm^3.
% 40 mm^3 is 3 voxels in all directions in DARTEL space, or around the highest obtainable ASL effective resolution
% (3x3x4)

end
