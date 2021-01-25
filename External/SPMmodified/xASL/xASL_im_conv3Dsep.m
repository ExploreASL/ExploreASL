function imConv = xASL_im_conv3Dsep(im, kX, varargin)
% xASL_im_conv3Dsep 3D separable convolution with a supplied kernel
% It converts the results to double
% Returned is the convoluted image
% The wrapper makes sure that kX are Nx1 format, removes nan, and removes excessive zeros at the ends
%
% FORMAT: [imConv] = xASL_mex_conv3Dsep(im,kX,[kY,kZ])
%
% INPUT:
% im - 3D image, double
% kX, kY, kZ - are the 1D kernels for the convolution
% 	- kY and kZ are optional
%      - if 0 or [] is supplied, then convolution in this dimension will be skipped
%
% OUTPUT:
% imConv - convolved image
%
% DESCRIPTION: 3D separable convolution with a supplied kernel
% It converts the results to double
% Returned is the convoluted image
% The wrapper makes sure that kX are Nx1 format, removes nan, and removes
% excessive zeros at the ends.
% 
% EXAMPLE: n/a
% __________________________________
% Copyright 2015-2021 ExploreASL

% Get the kY and kZ from the varargin
if nargin < 1
	kY = [];
end

if nargin < 2
	kZ = [];
end

if nargin > 0
	kY = varargin{1};
end

if nargin > 1
	kZ = varargin{2};
end

% Remove the nans
kX(isnan(kX)) = 0;
kY(isnan(kY)) = 0;
kZ(isnan(kZ)) = 0;

kX = kX(:);
kY = kY(:);
kZ = kZ(:);

% Length of the kernel must be odd
if (~isempty(kX)) && (~mod(length(kX),2))
	error('xASL_im_conv3Dsep: Only accepting odd-sized kernels');
end

if (~isempty(kY)) && (~mod(length(kY),2))
	error('xASL_im_conv3Dsep: Only accepting odd-sized kernels');
end

if (~isempty(kZ)) && (~mod(length(kZ),2))
	error('xASL_im_conv3Dsep: Only accepting odd-sized kernels');
end

% Remove excessive zeros at the end, while keeping the kernel center in center
% And normalize to 1
if ~isempty(kX)
	winSize = (length(kX)-1)/2;
	% Find the closest non-zero from the tail
	minNonzero = min([find(kX(1:winSize)~=0);find(kX(end:-1:(end-winSize+1))~=0)]);
	% Remove the zeros if it makes sense
	if (~isempty(minNonzero)) && (minNonzero>1)
		kX = kX(minNonzero:(end-minNonzero+1));
	end
	
	kX = kX/sum(kX(:));
end

if ~isempty(kY)
	winSize = (length(kY)-1)/2;
	minNonzero = min([find(kY(1:winSize)~=0);find(kY(end:-1:(end-winSize+1))~=0)]);
	if (~isempty(minNonzero)) && (minNonzero>1)
		kY = kY(minNonzero:(end-minNonzero+1));
	end
	
	kY = kY/sum(kY(:));
end

if ~isempty(kZ)
	winSize = (length(kZ)-1)/2;
	minNonzero = min([find(kZ(1:winSize)~=0);find(kZ(end:-1:(end-winSize+1))~=0)]);
	if (~isempty(minNonzero)) && (minNonzero>1)
		kZ = kZ(minNonzero:(end-minNonzero+1));
	end
	
	kZ = kZ/sum(kZ(:));
end

% Run the mex file
imConv = xASL_mex_conv3Dsep(im,kX,kY,kZ);
return;