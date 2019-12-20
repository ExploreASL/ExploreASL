function [imPVEC,imCBFrec,imResidual] = xASL_im_PVCkernel(imCBF, imPV,kernel,mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PVEc correction of ASL data using prior GM-,WM-partial volume maps.
% Follows the principles of the PVEc algorithm by I. Asllani (MRM, 2008).
% Free for research use without guarantee. If used in a study or
% publication. Please, acknowledge the author.
% Created by Jan Petr, j.petr@hzdr.de
% Version 0.8, 2018-06
%
% Please cite:
% -Asllani I, Borogovac A, Brown TR. Regression algorithm correcting for 
%  partial volume effects in arterial spin labeling MRI. Magnetic Resonance 
%  in Medicine. 2008 Dec 1;60(6):1362-71.
% -Petr J, Mutsaerts HJ, De Vita E, Steketee RM, Smits M, Nederveen AJ, 
%  Hofheinz F, van den Hoff J, Asllani I. Effects of systematic partial 
%  volume errors on the estimation of gray matter cerebral blood flow with 
%  arterial spin labeling MRI. MAGMA 2018. DOI:10.1007/s10334-018-0691-y
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   imCBF - uncorrected CBF image imCBF(X,Y,Z)
%   imPV  - PV-maps imPV(X,Y,Z,2) - WM/GM does not matter, you get the same
%           order in output as it is defined for input
%   kernel - 3D window size [I,J,K]
%   mode   - 'asllani' - standard implementation, flat kernel with size
%             [I,J,K], all the numbers must be odd.
%            'gauss' - Gaussian kernel is used with [I,J,K] is FWHM in
%             voxels, the window size is adjusted automatically to match
%             the Gaussians
%
% OUTPUT: 
%   imPVEC - corrected CBF maps for both tissues - imPVEC(X,Y,Z,2)
%   imCBFrec - (X,Y,Z) - reconstruction of CBF from the corrected values
%                        and PV maps
%   imResidual(X,Y,Z) - difference between the reconstructed and original
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check the input parameters
switch (mode)
	case 'asllani'
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
		error('xASL_im_PVCkernel: Only gauss and asllani modes are supported.');
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

imTotPV = sum(imPV,4);

% Creates a mask of the region covered by tissue
imMask = (imTotPV > opt.pvTh);

imCBF = imCBF.*imMask;
imPV = imPV.*repmat(imMask,[1 1 1 2]);

% Setting of the parameters for the lsqlin function
optLsqlin = optimset('Display','off','TolX',0.01,'TolFun',0.01,'FunValCheck','off','MaxIter',30);

% Specify the vectors for defining the kernel
xVec = (-winWidth(1)):winWidth(1);
yVec = (-winWidth(2)):winWidth(2);
zVec = (-winWidth(3)):winWidth(3);

% Creates the empty images for the tissue specific magnetization
imPVEC = zeros(size(imPV));

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
			imPVEC(x,y,z,1:2) = [0,0];
		
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
							
							imPVEC(x,y,z,1:2) = res;
							
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
								imPVEC(x,y,z,indMin) = opt.resMin;
								
								% Calculate the result for the other value
								imPVEC(x,y,z,indMax) = sum(ASLMatUp.*PMat(:,indMax))/sum(PMat(:,indMax).^2);
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
							
							imPVEC(x,y,z,ind) = sum(ASLMat.*PMat(:,ind))/sum(PMat(:,ind).^2);
						end
					end
				end
			end
		end
	end
end
%toc

% Reconstructed CBF based on CBF-PVEC and PV maps
imCBFrec = sum(imPV.*imPVEC,4);
        
% The residual error
imResidual = imMask.*(abs(imCBFrec-imCBF));

return;
