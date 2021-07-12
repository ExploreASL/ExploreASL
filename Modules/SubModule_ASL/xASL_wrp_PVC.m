function xASL_wrp_PVC(x)
% xASL_wrp_PVC Submodule of ExploreASL ASL Module, to perform PV correction in ASL native space
%
% FORMAT:  xASL_wrp_PVC(x)
%
% INPUT:
%   x - structure containing fields with all information required to run this submodule (REQUIRED)
%       x.modules.asl.PVCNativeSpaceKernel - Window size for the ASL native space PV correction. Equal weighting 
%                                         of all voxels within the kernel is assumed. 3D kernel can be used, 
%                                         but any of the dimension can be also set to 1. Only odd number of 
%                                         voxels can be used in each dimension (e.g. [3 7 5] not [2 3 1]).
%                                         (OPTIONAL, 
%                                         DEFAULT = [5 5 1] for bPVCGaussianMM==0,
%                                         DEFAULT = [10 10 4] for bPVCGaussianMM==1).
%       x.bPVCGaussianMM - PV-correction with a Gaussian instead of square kernel. It uses Gaussian 
%                                   weighting of the PV kernel instead of equal weights as per Asllani's 
%                                   original method. Unlike with the square kernel when the size is defined in 
%                                   voxels, here the FWHM of the Gaussian in mm is defined in each dimension. 
%                                   The advantage is twofold - continuous values can be added and a single 
%                                   value can be entered which is valid for datasets with different 
%                                   voxel-sizes without having a kernel of different effective size.
%                                   (OPTIONAL, DEFAULT = 0)
%                                    1 = enabled, use Gaussian kernel with FWHM in mm given in PVCNativeSpaceKernel
%                                    0 = disabled, use flat kernel with voxels given in PVCNativeSpaceKernel
%
% OUTPUT FILES: NIfTI containing PV-corrected CBF maps in ASL resolution in ASL native space
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This submodule performs partial volume correction (PVC) in native ASL space. It runs the Asllani's method
%              for partial volume correction by linear regression. It has two main extensions - first it uses a 3D kernel.
%              Second, it can use a Gaussian weights instead of the default flat kernel.
%              
%              0. Admin and checking values and files
%              1. Getting the resolution and preparing parameters
%              2. Running PV-correction
%              3. Saving files and cleaning
%
%
% EXAMPLE: xASL_wrp_PVC(x);
% __________________________________
% Copyright (C) 2015-2021 ExploreASL

% -----------------------------------------------------------------------------------------------------------------------------------------------------
%% 0. Admin and checking values and files
if ~isfield(x,'bPVCGaussianMM') || isempty(x.bPVCGaussianMM)
	x.bPVCGaussianMM = 0;
end

% If the kernel is non-existent or empty then initialize it
if ~isfield(x.modules.asl,'PVCNativeSpaceKernel') || isempty(x.modules.asl.PVCNativeSpaceKernel)
	if x.bPVCGaussianMM
		x.modules.asl.PVCNativeSpaceKernel = [10 10 4];
	else
		x.modules.asl.PVCNativeSpaceKernel = [5 5 1];
	end
end

% If the size of the kernel was under 3, then add the remaining dimensions from the default
dimKernel = length(x.modules.asl.PVCNativeSpaceKernel);
if dimKernel < 3
	if x.bPVCGaussianMM
		defaultKernel = [10 10 4];
	else
		defaultKernel = [5 5 1];
	end
	x.modules.asl.PVCNativeSpaceKernel(dimKernel:3) = defaultKernel(dimKernel:3);
end

% -----------------------------------------------------------------------------------------------------------------------------------------------------
%% 1. Getting the resolution and preparing parameters
% Load the CBF image
if ~xASL_exist(x.P.Path_CBF,'file')
	error(['CBF file ' x.P.Path_CBF ' does not exist. Skipping PVC']);
end

imCBF = xASL_io_Nifti2Im(x.P.Path_CBF);

% Load and concatenate the PV maps
if ~xASL_exist(x.P.Path_PVwm,'file')
	error(['PV-WM file ' x.P.Path_PVwm ' does not exist. Skipping PVC']);
end
if ~xASL_exist(x.P.Path_PVgm,'file')
	error(['PV-GM file ' x.P.Path_PVgm ' does not exist. Skipping PVC']);
end

imGM = xASL_io_Nifti2Im(x.P.Path_PVgm);
imWM = xASL_io_Nifti2Im(x.P.Path_PVwm);

% Check if size is compatible
if ~isequal(size(imGM),size(imWM)) || ~isequal(size(imGM),size(imCBF))
	error('Size of the GM, WM, and CBF files have to match.');
end

imPV = imGM;
imPV(:,:,:,2) = imWM;

if x.bPVCGaussianMM
	% Prepare the kernel size. Convert the FWHM from voxels to MM
	voxelSize = xASL_io_ReadNifti(x.P.Path_CBF);
	voxelSize = [norm(voxelSize.mat(:,1)), norm(voxelSize.mat(:,2)), norm(voxelSize.mat(:,3))];
	
	kernelPVC = x.modules.asl.PVCNativeSpaceKernel./voxelSize;
else
	kernelPVC = x.modules.asl.PVCNativeSpaceKernel;
end

% -----------------------------------------------------------------------------------------------------------------------------------------------------
%% 2. Running PV-correction

if x.bPVCGaussianMM
	[imPVC,~,~] = xASL_im_PVCkernel(imCBF, imPV, kernelPVC, 'gauss');
else
	[imPVC,~,~] = xASL_im_PVCkernel(imCBF, imPV, kernelPVC, 'flat');
end

% -----------------------------------------------------------------------------------------------------------------------------------------------------
%% 3. Saving files and cleaning
% Remove areas of GM-CBF and WM-CBF with respective GM-PV and WM-PV under 3 percent
imPVC(:,:,:,1) = imPVC(:,:,:,1).*(imGM>0.1);
imPVC(:,:,:,2) = imPVC(:,:,:,2).*(imWM>0.1);
imPVC(:,:,:,1) = imPVC(:,:,:,1).*((imGM+imWM)>0.3);
imPVC(:,:,:,2) = imPVC(:,:,:,2).*((imGM+imWM)>0.3);

% Save the GM-CBF and WM-CBF files to the individual ASL directory
xASL_io_SaveNifti(x.P.Path_CBF,x.P.Path_CBFgm,imPVC(:,:,:,1));
xASL_io_SaveNifti(x.P.Path_CBF,x.P.Path_CBFwm,imPVC(:,:,:,2));

end
