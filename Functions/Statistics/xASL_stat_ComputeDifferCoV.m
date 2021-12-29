function diffCoV = xASL_stat_ComputeDifferCoV(imCBF,imMask,bPVC,imGM,imWM,b3D)
% Calculates spatial diff coefficient of variation (sCoV) in the image with optional partial volume correction.
%
% FORMAT: diffCoV = xASL_stat_ComputeDifferCoV(imCBF)
%         diffCoV = xASL_stat_ComputeDifferCoV(imCBF,imMask)
%         diffCoV = xASL_stat_ComputeDifferCoV(imCBF,imMask,bPVC,imGM,imWM)
%         diffCoV = xASL_stat_ComputeDifferCoV(imCBF,imMask,bPVC,imGM,imWM,b3D)
%
% INPUT:
%   imCBF  - input CBF volume
%   imMask - mask for the calculation (DEFAULT finite part of imCBF)
%   bPVC   - perform PV-correction (DEFAULT 0)
%            0 - don't do partial volume correction
%            2 - partial volume correction by using pseudoCov calculated from imGM, imWM
%   imGM   - GM partial volume map with the same size as imCBF, mandatory for bPVC==2
%   imWM   - WM partial volume map with the same size as imCBF, mandatory for bPVC==2
%   b3D    - calculate 2D wise or 3D wise (DEFAULT 0)
% OUTPUT:
%   sCov   - calculated spatial coefficient of variation
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: It calculates the spatial DiffCoV value on finite part of imCBF. Optionally a mask IMMASK is provide, 
%              and PVC is done for bPVC==2 using imGM and imWM masks and constructing
%              pseudoCoV from pseudoCBF image. For bPVC~=2, imGM and imWM are ignored. It is calculated in 2D or assuming also 3D edges based on B3D.
%              Calculate derivate spatial CoV, by summing up differences in CBF between neighbors.
%              The derivative uses Sobels filter.
%
% EXAMPLE: diffCoV = xASL_stat_ComputeDifferCoV(imCBF)
%          diffCoV = xASL_stat_ComputeDifferCoV(imCBF,imMask)
%          diffCoV = xASL_stat_ComputeDifferCoV(imCBF,[],0)
%          diffCoV = xASL_stat_ComputeDifferCoV(imCBF,imMask,2,imGM,imWM)
%          diffCoV = xASL_stat_ComputeDifferCoV(imCBF,imMask,2,imGM,imWM,1)
%          diffCoV = xASL_stat_ComputeDifferCoV(imCBF,[],0,[],[],1)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright © 2015-2019 ExploreASL
%
% 2017-00-00 JP

% Admin
if nargin < 2 || isempty(imMask)
	imMask = ones(size(imCBF));
end

if nargin < 3 || isempty(bPVC)
	bPVC = 0;
end

if nargin < 4
	imGM = [];
end

if nargin < 5
	imWM = [];
end

if nargin < 6 || isempty(b3D)
	b3D = 0;
end

if bPVC < 2
	if ~isempty(imGM) || ~isempty(imWM)
		error('xASL_stat_ComputeDifferCoV: imGM and imWM should only be provided for bPVC==2.');
	end
end

if bPVC == 2
	% If running PVC, then need imGM and imWM of the same size as imCBF
	if (~isequal(size(imCBF),size(imGM))) || (~isequal(size(imCBF),size(imWM)))
		error('xASL_stat_ComputeDifferCoV: When running PVC, need imGM and imWM of the same size as imCBF');
	end
end

% Initialize the filters
sobelFilter     = [-1 -2 -1;0 0 0; 1 2 1];
averageFilter   = [1 1 1; 1 1 1; 1 1 1];
tempEdge        = zeros(size(imCBF),'single');
tempAverage     = zeros(size(imCBF),'single');

% Calculate the pseudoCoV as well for partial volume normalization
if bPVC == 2
	pseudoCoV    = zeros(size(imCBF),'single');
	imPseudoCBF  = imGM + 0.3.*imWM;
end
    
% Do for all slices and sum for both in-plane directions
for sliceIdx = 1:size(imCBF,3)
	tempEdge(:,:,sliceIdx)      = (conv2(imCBF(:,:,sliceIdx),sobelFilter,'same')/8).^2 + (conv2(imCBF(:,:,sliceIdx),sobelFilter','same')/8).^2;
	tempAverage(:,:,sliceIdx)   = conv2(imCBF(:,:,sliceIdx),averageFilter,'same');
        
	% Calculates the average over the same mask to make sure that only reasonable values
	% are taken and no nans from outside are used
        
	if  bPVC==2
		pseudoCoV(:,:,sliceIdx) = (conv2(imPseudoCBF(:,:,sliceIdx),sobelFilter,'same')/8).^2 + (conv2(imPseudoCBF(:,:,sliceIdx),sobelFilter','same')/8).^2;
		tempAverage(:,:,sliceIdx)   = tempAverage(:,:,sliceIdx) + conv2(imPseudoCBF(:,:,sliceIdx),averageFilter,'same');
	end
end
    
if b3D
	% Do the same for the coronal and sagittal slices
	for sliceIdx = 1:size(imCBF,2)
		sliceData = squeeze(imCBF(:,sliceIdx,:));
		tempEdge(:,sliceIdx,:)      = tempEdge(:,sliceIdx,:) + (conv2(sliceData,sobelFilter,'same')/8).^2 + (conv2(sliceData,sobelFilter','same')/8).^2;
		tempAverage(:,sliceIdx,:)   = tempAverage(:,sliceIdx,:) + conv2(sliceData,averageFilter,'same');
            
		if bPVC == 2
			sliceData = squeeze(imPseudoCBF(:,sliceIdx,:));
			pseudoCoV(:,sliceIdx,:) = pseudoCoV(:,sliceIdx,:) + (conv2(sliceData,sobelFilter,'same')/8).^2 + (conv2(sliceData,sobelFilter','same')/8).^2;
			tempAverage(:,sliceIdx,:)   = tempAverage(:,sliceIdx,:) + conv2(sliceData,averageFilter,'same');
		end
	end
        
	for sliceIdx = 1:size(imCBF,1)
		sliceData = squeeze(imCBF(sliceIdx,:,:));
		tempEdge(sliceIdx,:,:)      = tempEdge(sliceIdx,:,:) + (conv2(sliceData,sobelFilter,'same')/8).^2 + (conv2(sliceData,sobelFilter','same')/8).^2;
		tempAverage(sliceIdx,:,:)   = tempAverage(sliceIdx,:,:) + conv2(sliceData,averageFilter,'same');
		
		if bPVC == 2
			sliceData = squeeze(imPseudoCBF(sliceIdx,:,:));
			pseudoCoV(sliceIdx,:,:) = pseudoCoV(sliceIdx,:,:) + (conv2(sliceData,sobelFilter,'same')/8).^2 + (conv2(sliceData,sobelFilter','same')/8).^2;
			tempAverage(sliceIdx,:,:)   = tempAverage(sliceIdx,:,:) + conv2(sliceData,averageFilter,'same');
		end
	end
end
    
tempEdge = sqrt(tempEdge);
pseudoCoV = sqrt(pseudoCoV);
    
imMask                  = imMask & isfinite(tempAverage);
    
% Calculate a mean of all data and of the derivatives
diffCoV         = mean(tempEdge(imMask)) / mean(imCBF(imMask));

if  bPVC==2
	pseudoCoV       = mean(pseudoCoV(imMask)) / mean(imPseudoCBF(imMask));
	diffCoV     = diffCoV./pseudoCoV;
end
    
end

