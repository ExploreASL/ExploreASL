function pathOut = xASL_im_PreSmooth(pathRef,pathSrc,pathSmo,resRef,resSrc,srcAffinePath, bInvAffine)
% xASL_im_PreSmooth takes the source image and pre-smooths to the lower resolution of the reference image,
% so that it can be easily downsampled using a normal trilinear interpolation without loosing information.
%
% FORMAT: pathOut = xASL_im_PreSmooth(pathRef,pathSrc[,pathSmo,resRef,resSrc,srcAffinePath, bInvAffine])
%
% INPUT:
%   pathRef        The path to the reference image with the final resolution/voxels. (REQUIRED)
%   pathSrc        The path to the source image to be smoothed to fit the reference image. (REQUIRED)
%   pathSmo        The path for saving the result. (OPTIONAL. DEFAULT adds a prefix 'r' to the source image)
%   resRef         Reference image resolution in mm. 3x1 vector. (OPTIONAL. DEFAULT matches the voxel size)
%   resSrc         Source image resolution in mm. 3x1 vector. (OPTIONAL. DEFAULT matches the voxel size)
%   srcAffinePath  The path to the '*_sn.mat' that describes the Affine transformation matrix related
%                  to the source image and describing the src->ref transformation.
%   bInvAffine     Do the inverse of the affine transform (OPTIONAL, default = FALSE)
%                  if the affine goes from src->ref, then = FALSE. For ref->src, we need to use the inversion = TRUE
%
% OUTPUT:
%   pathOut        The path to the output file.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: It assumes that the FWHM is equal to voxel size, unless the real resolution is given.
%              Then takes into account the voxel sizes and orientation difference between the volumes, but
%              performs the smoothing according to the given real resolution (it is possible to supply the
%              resolution for just one image) - this should be helpful primarily when the
%              It creates a Gaussian covariance matrix for the reference, transforms it according to
%              the difference between the two images a produces the Gaussian covariance matrix
%              of the pre-smoothing in the source space. Then it performs the smoothing.
% EXAMPLE:
%     pathOut = xASL_im_PreSmooth('/home/tmp/CBF.nii','/home/tmp/T1.nii')
%     pathOut = xASL_im_PreSmooth('/home/tmp/CBF.nii','/home/tmp/T1.nii',[],[3 3 12],[])
%     pathOut = xASL_im_PreSmooth('/home/tmp/CBF.nii','/home/tmp/T1.nii',[],[3 3 12],[],'/home/tmp/CBF_sn.mat',1)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% REFERENCES:
% Pre-smoothing of the source image before interpolation according to:
% Cardoso M.J., Modat M., Vercauteren T., Ourselin S. (2015) Scale Factor Point Spread Function Matching:
% Beyond Aliasing in Image Resampling. In: Navab N., Hornegger J., Wells W., Frangi A. (eds) Medical Image
% Computing and Computer-Assisted Intervention -- MICCAI 2015. MICCAI 2015. Lecture Notes in Computer Science,
% vol 9350. Springer, Cham
% __________________________________
% Copyright © 2015-2019 ExploreASL
%
% 2019-06-14 JP

% Admin
if nargin < 2 || isempty(pathRef) || isempty(pathSrc)
	error('xASL_im_PreSmooth: Need to specify the reference and source images.');
end

% Create the output by adding a prefix 'r'
if nargin < 3 || isempty(pathSmo)
	[Fpath, Ffile, Fext] = xASL_fileparts(pathSrc);
	pathOut = fullfile(Fpath,['r' Ffile Fext]);
else
	pathOut = pathSmo;
end

% Loads the input images
imRef = xASL_io_ReadNifti(pathRef);
imSrc = xASL_io_ReadNifti(pathSrc);

% Define the reference and source resolution as FWHM in voxels
if nargin < 4 || isempty(resRef)
	% By default equal to the voxel size
	resVoxRef = [1 1 1];
else
	% Obtain the voxel size and divide the resolution by the voxel size
	voxRef = [norm(imRef.mat(:,1)), norm(imRef.mat(:,2)), norm(imRef.mat(:,3))];
	resVoxRef = resRef./voxRef;
end

if nargin < 5 || isempty(resSrc)
	% By default equal to the voxel size
	resVoxSrc = [1 1 1];
else
	% Obtain the voxel size and divide the resolution by the voxel size
	voxSrc = [norm(imSrc.mat(:,1)), norm(imSrc.mat(:,2)), norm(imSrc.mat(:,3))];
	resVoxSrc = resSrc./voxSrc;
end

if min(resVoxSrc==resVoxRef)
    xASL_Copy(pathSrc, pathOut, 1);
    return; % skip pre-smoothing
end


if nargin < 6 || isempty(srcAffinePath)
	srcAffinePath = [];
elseif iscell(srcAffinePath)
	srcAffinePath = srcAffinePath{1};
end

if nargin < 7 || isempty(bInvAffine)
	bInvAffine = false;
end

% Checks if the sn-mat exists and then loads it
if ~isempty(srcAffinePath) && exist(srcAffinePath,'file')
	SnMat = load(srcAffinePath);
else
	SnMat = [];
end

% Check if the affine transform was loaded and then modify the source matrix accordingly
if isempty(SnMat)
	srcMat = imSrc.mat;
else
	if bInvAffine
		srcMat = SnMat.VF.mat*SnMat.Affine;
	else
		srcMat = SnMat.VG.mat/SnMat.Affine;%SnMat.VG.mat*inv(SnMat.Affine);
	end
	% There is an external Mat file provided or potentially different internal MAT (different from the Affine trans definition) - take this into account
	if bInvAffine
		% Doing the inverse
		srcMat = (srcMat/SnMat.VG.mat)*imSrc.mat;%mat*inv(SnMat.VG.mat)*locMat;
	else
		srcMat = (srcMat/SnMat.VF.mat)*imSrc.mat;%mat*inv(SnMat.VF.mat)*locMat;
	end
end

% Obtain the transformation matrix from the Reference to the Source space
%transMat = inv(imSrc.mat)*imRef.mat;
transMat = srcMat\imRef.mat;
transMat = transMat(1:3,1:3);

% The covariance matrix in both spaces is simply a diagonal matrix of the resolution in voxels
% Multiply by the transformation matrix to bring the covariance matrix of the reference space to the source
% space and subtract the covariance matrix of the source to obtain the covariance matrix of the filter
covMat = (transMat*diag(resVoxRef.^2)*(transMat') - diag(resVoxSrc.^2));

% covMat is in FWHM^2 we need to move it to var by dividing by (2*sqrt(2*log(2)))^2
covMat = covMat/((2*sqrt(2*log(2)))^2);

% The inverse of the covariance matrix is needed to calculate the 3D Gaussian
%covMatInv = inv(covMat);

% Obtain the maximal sigma. We cut the window at 2.7 sigma.
wSize = ceil(sqrt(max(abs(covMat),[],2))*2.7);

% Prepare the pixel grid for calculating of the kernel
[pX,pY,pZ] = ndgrid(-wSize(1):wSize(1),-wSize(2):wSize(2),-wSize(3):wSize(3));
pVec = [pX(:),pY(:),pZ(:)]';

% Prepare the empty kernel
kFil = zeros(size(pVec,2),1);
for ii=1:size(pVec,2)
	%kFil(ii) = pVec(:,ii)'*covMatInv*pVec(:,ii);
	kFil(ii) = (pVec(:,ii)'/covMat)*pVec(:,ii);
end
% Reshape back to normal size
kFil = reshape(kFil,size(pX));

% Calculate the PDF
kFil = exp(-kFil/2);

% Normalize the PDF
kFil = kFil/sum(kFil(:));

% Apply the smoothing filter on the source image
imSmo = zeros(size(imSrc.dat));
for ii = 1:size(imSrc.dat,4)
	imSmo(:,:,:,ii) = convn(imSrc.dat(:,:,:,ii),kFil,'same');
end

% Save the smoothed image
xASL_io_SaveNifti(pathSrc, pathOut, imSmo, [], 0);

end
