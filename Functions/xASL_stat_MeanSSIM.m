function mssim=xASL_stat_MeanSSIM(imRef,imSrc,dynRange)
% Calculates the mean structural similarity index 
%
% FORMAT: mssim=xASL_stat_MeanSSIM(imRef,imSrc[,dynRange])
%
% INPUT:
%   imRef    - input reference image
%   imSrc    - input source image
%   dynRange - dynamic range of the images used for the meanSSIM calculation
%            - by default takes the maximum from the reference image
% OUTPUT:
%   mmsim    - mean SSIM
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Calculates the similarity index according to Want et al. 
%
% EXAMPLE: xASL_stat_MeanSSIM(imRef, imSrc, 255)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCE: Z. Wang, A. C. Bovik, H. R. Sheikh and E. P. Simoncelli, Image quality assessment: From error 
%            visibility to structural similarity, IEEE Transactions on Image Processing, 13(4):600-612, 2004.
% __________________________________
% Copyright (C) 2015-2019 ExploreASL
%
% 2019-04-24 JP, HM

% Initialization
if nargin < 2
	error('xASL_stat_MeanSSIM: Needs to input images.');
end

if ~isequal(size(imRef),size(imSrc))
	error('xASL_stat_MeanSSIM: Input images need to have the same size.');
end

if nargin < 3 || isempty(dynRange)
	% By default, set the dynamic range to maximal values
	dynRange = max(abs(min(imRef(:))),abs(max(imRef(:))));
end

% Constants set according to the original article
K1 = 0.01;
K2 = 0.03;

C1 = (K1*dynRange).^2;                  
C2 = (K2*dynRange).^2;

imRef = double(imRef);
imSrc = double(imSrc);

% Mean over the local window and its square
mu1         = xASL_im_ndnanfilter(imRef,'gauss',[4 4 4]);
mu2         = xASL_im_ndnanfilter(imSrc,'gauss',[4 4 4]);
mu1_2       = mu1.^2;
mu2_2       = mu2.^2;
mu1_mu2     = mu1.*mu2;

% Squared standard deviation computed as mean of square - square of mean
sigma1_2    = xASL_im_ndnanfilter(imRef.^2,'gauss',[4 4 4]) - mu1_2;
sigma2_2    = xASL_im_ndnanfilter(imSrc.^2,'gauss',[4 4 4]) - mu2_2;
sigma12     = xASL_im_ndnanfilter(imRef.*imSrc,'gauss',[4 4 4]) - mu1_mu2;

% Calculate the similarity map and the mean SSIM
ssim_map    =  ((2*mu1_mu2 + C1).*(2*sigma12 + C2)) ./ ((mu1_2 + mu2_2 + C1).*(sigma1_2 + sigma2_2 + C2));
mssim       = 100*xASL_stat_MeanNan(ssim_map(:)); % in percentages

end
