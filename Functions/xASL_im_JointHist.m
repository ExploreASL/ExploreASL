function imHist = xASL_im_JointHist(imA,imB,imMask,minA,maxA,minB,maxB,nBins)
% xASL_im_JointHist calculates a joint histogram of two images across a binary mask
%
% FORMAT: imHist = xASL_im_JointHist(imA,imB[,imMask,minA,maxA,minB,maxB,nBins])
%
% INPUT:
%   imA    - First input image (REQUIRED).
%   imB    - Second input image, needs to have the same dimensions (REQUIRED).
%   imMask - Binary mask to calculate the joint histogram (OPTIONAL, DEFAULT whole image).
%   minA   - Minimal value to be counted for image A (OPTIONAL, DEFAULT imA minimum).
%   maxA   - Maximal value to be counted for image A (OPTIONAL, DEFAULT imA maximum).
%   minB   - Minimal value to be counted for image B (OPTIONAL, DEFAULT imB minimum).
%   minB   - Maximal value to be counted for image B (OPTIONAL, DEFAULT imB maximum).
%   nBins  - Number of bins (OPTIONAL, DEFAULT 20).
%
% OUTPUT:
%   imHist - Resulting joint histogram of size nBins x nBins.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: It calculates a joint histogram of two images of any dimensions over a mask of the same size.
%              The boundaries and number of bins can either be given or min and max values are used. Values
%              outside of the bins are counted to the first/last bin.
% EXAMPLE: 
%     imHist = xASL_im_JointHist(imA,imB)
%     imHist = xASL_im_JointHist(imA,imB,[],0,200,-50,[])
%     imHist = xASL_im_JointHist(imA,imB,imMask,[],200,0,[],50)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright © 2015-2019 ExploreASL
%
% 2019-06-26 JP

% Admin
if nargin < 2
	error('xASL_im_JointHist: Needs at least two images to calculate histogram from.');
end

if ~isequal(size(imA),size(imB))
	error('xASL_im_JointHist: The two input images need to have the same size.');
end

if nargin < 3 || isempty(imMask)
	imMask = ones(size(imA));
else
	if ~isequal(size(imA),size(imMask))
		error('xASL_im_JointHist: The input images and the mask need to have the same size.');
	end
end

if nargin < 4 || isempty(minA)
	minA = min(imA(:));
end

if nargin < 5 || isempty(maxA)
	maxA = min(imA(:));
end

if nargin < 6 || isempty(minB)
	minB = min(imB(:));
end

if nargin < 7 || isempty(maxB)
	maxB = min(imB(:));
end

if nargin < 8 || isempty(nBins)
	nBins = 20;
end

imHist = xASL_mex_JointHist(imA,imB,imMask,minA,maxA,minB,maxB,nBins);

end