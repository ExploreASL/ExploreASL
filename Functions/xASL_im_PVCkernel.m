function [imPVC,imCBFrec,imResidual] = xASL_im_PVCkernel(imCBF, imPV,kernel,mode)
%xASL_im_PVCkernel Partial volume correction (PVC) of ASL data using GM-,WM-partial volume maps using linear regression
%
% FORMAT: [imPVC,imCBFrec,imResidual] = xASL_im_PVCkernel(imCBF, imPV [,kernel,mode])
%
% INPUT:
%   imCBF  - uncorrected CBF image imCBF(NX,NY,NZ) (REQUIRED)
%   imPV   - PV-maps imPV(NX,NY,NZ,2) - WM/GM order does not matter, you get the same
%            order in output as it is defined for input (REQUIRED)
%   kernel - 3D window size [NI,NJ,NK] (OPTIONAL, DEFAULT = [5 5 1])
%   mode   - Type of the kernel used - flat or Gaussian (OPTIONAL, DEFAULT = 'flat')
%            'flat' - standard implementation, flat kernel with size
%             [NI,NJ,NK], all the numbers must be odd.
%            'gauss' - Gaussian kernel is used with [NI,NJ,NK] is FWHM in
%             voxels, the window size is adjusted automatically to match
%             the Gaussians
%
% OUTPUT: 
%   imPVC - corrected CBF maps for both tissues - imPVC(NX,NY,NZ,2)
%   imCBFrec - (NX,NY,NZ) - reconstruction of CBF from the corrected values
%                        and PV maps
%   imResidual(NX,NY,NZ) - difference between the reconstructed and original
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Partial volume correction (PVC) of ASL data using prior GM-,WM-partial volume maps.
%               Follows the principles of the PVC algorithm by I. Asllani (MRM, 2008).
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      [imPVC,imCBFrec,imResidual] = xASL_im_PVCkernel(imCBF, imPV,[5 5 1],'flat')
%
% REFERENCES:
% -Asllani I, Borogovac A, Brown TR. Regression algorithm correcting for 
%  partial volume effects in arterial spin labeling MRI. Magnetic Resonance 
%  in Medicine. 2008 Dec 1;60(6):1362-71.
% -Petr J, Mutsaerts HJ, De Vita E, Steketee RM, Smits M, Nederveen AJ, 
%  Hofheinz F, van den Hoff J, Asllani I. Effects of systematic partial 
%  volume errors on the estimation of gray matter cerebral blood flow with 
%  arterial spin labeling MRI. MAGMA 2018. DOI:10.1007/s10334-018-0691-y
% __________________________________
% Copyright 2015-2021 ExploreASL

%% 0. Admin
if nargin<2
	error('CBF image and PV-maps have to be given on the input');
end

if nargin<3 || isempty(kernel)
	kernel = [5 5 1];
end

if nargin<4 || isempty(mode)
	mode = 'flat';
end

%% 1. Prepare the parameters
% Check the input parameters
switch (lower(mode))
	case 'flat'
		% The kernel size has to be odd in all dimensions
		if (sum(mod(kernel,2))<3)
			error('xASL_im_PVCkernel: kernel size has to be odd in all dimensions.');
		end
		winWidth = (kernel-1)/2;
		doGauss = 0;
	case 'gauss'
		% Calculate the kernel width to fit the Gaussian
		winWidth = ceil(2.5*kernel/2.355);
		
		% Calculate the Gaussian in all dimensions
		for i=1:3
			imGauss{i} = exp( -((-winWidth(i):winWidth(i)).^2)./2./((kernel(i)/2.355)^2));
			imGauss{i} = imGauss{i}/max(imGauss{i});
		end
		
		% Create the Gaussian kernel in 2D
		imGaussKernel = imGauss{1}'*imGauss{2};
		
		% Create the 3D kernel
		imGaussKernel = repmat(imGaussKernel,[1 1 length(imGauss{3})]);
		imGaussZ = repmat(imGauss{3},[1 1 1]);
		imGaussZ = shiftdim(imGaussZ,-1);
		imGaussKernel = imGaussKernel .* repmat(imGaussZ,[length(imGauss{1}),length(imGauss{2}),1]);
		
		doGauss = 1;
	otherwise
		error('xASL_im_PVCkernel: Only gauss and flat modes are supported.');
end

imPV(isnan(imPV)) = 0;
imCBF(isnan(imCBF)) = 0;

% Tissue below this threshold has not PVEC calculated, nor is considered
% for the neighboring voxel PVEC calculation.
opt.pvTh = 0.01;

% The sum of PV across the whole kernel has to be above this to all the
% calculations of PVEC
opt.pvTotTh = 0.1;

% There has to be at least 1 pixel with partial volume higher than 'opt.rankTh' to compute the inversion
opt.rankTh = 0.01;

if doGauss
	opt.pvTh = 0.002;
	opt.pvTotTh = 0.5;
	opt.rankTh = 0.002;
end

% Any value below this is 
opt.resMin = 0;

%% 2. Performs the PV correction
imTotPV = sum(imPV,4);

% Creates a mask of the region covered by tissue
imMask = (imTotPV > opt.pvTh);

imCBF = imCBF.*imMask;
imPV = imPV.*repmat(imMask,[1 1 1 2]);

% Specify the vectors for defining the kernel
xVec = (-winWidth(1)):winWidth(1);
yVec = (-winWidth(2)):winWidth(2);
zVec = (-winWidth(3)):winWidth(3);

% Creates the empty images for the tissue specific magnetization
imPVC = zeros(size(imPV));

%tic
% For each pixel on the mask do the calculation
for z=1:size(imCBF,3)
	% Since the FOV in Z-direction is limited, we always go through all the
	% voxels in the Z-direction, but we have to adapt the kernel not to go
	% out of the volume
	
	zVecLoc = zVec;
	zVecLoc = zVecLoc((z+zVecLoc)>0);
	zVecLoc = zVecLoc((z+zVecLoc)<=size(imCBF,3));
	numElLoc = length(xVec) * length(yVec) * length(zVecLoc);
	
	for y = (winWidth(2)+1):(size(imCBF,2)-winWidth(2))
		for x = (winWidth(1)+1):(size(imCBF,1)-winWidth(1))
			imPVC(x,y,z,1:2) = [0,0];
		
			if imMask(x,y,z)
				% Create the empty vectors
				PMat = zeros(numElLoc,2);
				
				if doGauss
					% Gaussian kernel
					kernelLoc = imGaussKernel(xVec + winWidth(1) + 1, yVec + winWidth(2) + 1, zVecLoc + winWidth(3) + 1);
					
					% Create the PV matrix
					PMat(:,1) = reshape(imPV(x + xVec, y + yVec, z + zVecLoc,1).*kernelLoc,[numElLoc,1]);
					PMat(:,2) = reshape(imPV(x + xVec, y + yVec, z + zVecLoc,2).*kernelLoc,[numElLoc,1]);
					
					% And CBF vector
					ASLMat = reshape(imCBF(x + xVec, y + yVec, z + zVecLoc,1).*kernelLoc,[numElLoc,1]);
				else 
					% Flat kernel
					% Create the PV matrix
					PMat(:,1) = reshape(imPV(x + xVec, y + yVec, z + zVecLoc,1),[numElLoc,1]);
					PMat(:,2) = reshape(imPV(x + xVec, y + yVec, z + zVecLoc,2),[numElLoc,1]);
					
					% And CBF vector
					ASLMat = reshape(imCBF(x + xVec, y + yVec, z + zVecLoc,1),[numElLoc,1]);
				end
				% Remove the parts with low pvGM+pvWM
				ind = (sum(PMat,2) > opt.pvTh) & (ASLMat~=0);
				PMat = PMat(ind,:);
				ASLMat = ASLMat(ind,:);
				
				% There has to be a minimum amount of tissue
				if (sum(PMat(:)) > opt.pvTotTh)
					
					% Calculate if there is voxels for both tissue
					isPV1 = (sum(PMat(:,1) > opt.rankTh) >= 1);
					isPV2 = (sum(PMat(:,2) > opt.rankTh) >= 1);
					
					if ( isPV1 && isPV2 )
						% Compute the Matrix rank of PMat
						[PMatu,PMats,PMatv] = svd(PMat,0);
						PMatds = diag(PMats);
						PMattol = max(size(PMat)) * eps(max(PMatds));
						PMatr = sum(PMatds > PMattol);
						
						% The matrix is not ill-conditioned
						if PMatr >= 2
							% Calculates the pseudo inversion and the
							% solution
							PMatds = 1./PMatds(:);
							%res = ((PMatv.*PMatds.')*(PMatu'))*ASLMat;
							res = (PMatv*diag(PMatds)*PMatu')*ASLMat;
							
							imPVC(x,y,z,1:2) = res;
							
							% In case that the result for one or both
							% tissue is below minimal allowed value, we set
							% the smaller to the minimum and recalculate
							
							% One of the values is smaller than the MIN
							if sum(res<opt.resMin)
								% Find out the index of the bigger one and
								% the smaller one
								[~,indMax] = max(res);
								[~,indMin] = min(res);
								% Set one value to the min and subtract it
								ASLMatUp = ASLMat - PMat(:,indMin)*opt.resMin;
								
								% Set the results for the smaller value to
								% the minimum
								imPVC(x,y,z,indMin) = opt.resMin;
								
								% Calculate the result for the other value
								imPVC(x,y,z,indMax) = sum(ASLMatUp.*PMat(:,indMax))/sum(PMat(:,indMax).^2);
							end
						end
					else
						% Not enough voxels of both tissue
						if (isPV1 || isPV2)
							% But one tissue has enough, so
							% calculate for the first of second tissue
							if (isPV1)
								ind = 1;
							else
								ind = 2;
							end
							
							imPVC(x,y,z,ind) = sum(ASLMat.*PMat(:,ind))/sum(PMat(:,ind).^2);
						end
					end
				end
			end
		end
	end
end
%toc

%% 3. Creates additional outputs
% Reconstructed CBF based on CBF-PVEC and PV maps
imCBFrec = sum(imPV.*imPVC,4);
        
% The residual error
imResidual = imMask.*(abs(imCBFrec-imCBF));

end

