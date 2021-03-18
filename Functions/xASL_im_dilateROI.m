function xASL_im_dilateROI(PathIn, PathOut, minVolume)
%xASL_im_dilateROI Dilates a binary ROI with a 3x3 sphere until reaching a defined minimal volume
%
% FORMAT:       xASL_im_dilateROI(PathIn, PathOut)
% 
% INPUT:        PathIn    - Path to input NIfTI image (REQUIRED)
%               PathOut   - Path to output NIfTI image (OPTIONAL, DEFAULT = PathIn)
%               minVolume - The ROI is dilated until reaching this volume in ml (OPTIONAL, DEFAULT 40ml)
% OUTPUT:       n/a
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  The function loads a binary image from PathIn and if smaller than the defined volume (40 mL by default) it 
%               dilates it with a 3x3 sphere element until a minimal volume is reached. When it is small enough, it is saved to PathOut.
%               40 mm^3 is equal to 3 voxels in all directions in DARTEL space, or around the highest obtainable ASL effective resolution (3x3x4 mm).
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      xASL_im_dilateROI('test.nii', [], 40)
% __________________________________
% Copyright 2015-2021 ExploreASL

%% Admin
if nargin<2 || isempty(PathOut)
    PathOut = PathIn;
end

if nargin<3 || isempty(minVolume)
	minVolume = 40;
end

%% Loads the image and converts to a binary mask
IM      = xASL_io_Nifti2Im(PathIn);
IM      = xASL_im_ConvertMap2Mask(IM);

%% Obtains the voxel volume
HDR     = xASL_io_ReadNifti(PathIn);
pixdim  = prod(HDR.hdr.pixdim(2:4));
IMvol   = sum(IM(:)).*pixdim;

%% Iteratively dilate until the minimal volume value is reached
while  IMvol<minVolume
       IM = xASL_im_DilateErodeFull(IM,'dilate',xASL_im_DilateErodeSphere(2));
       IMvol   = sum(IM(:)).*pixdim;
end

%% Saves the dilated ROI
xASL_io_SaveNifti(PathIn, PathOut, IM, [], false);

end
