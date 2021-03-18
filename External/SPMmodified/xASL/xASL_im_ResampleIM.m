function [imOut] = xASL_im_ResampleIM(imIn, matIn, matOut, dimOut, interpolationType)
%xASL_im_ResampleIM Resamples an image with provided input and output orientation matrices
%
% FORMAT: [imOut] = xASL_im_ResampleIM(imIn, matIn, matOut, dimOut[, interpolationType])
%
% INPUT:
%   imIn	          - input image for resampling - an image matrix (REQUIRED)
%   matIn	          - 4x4 orientation matrix of the input image (REQUIRED)
%   matOut	          - 4x4 orientation matrix of the output image (REQUIRED)
%   dimOut	          - 4D vector of output image dimensions (REQUIRED)
%   interpolationType - Type of interpolation used by interp3 - 'linear' 'cubic' 'or 'nearest'
%                       (OPTIONAL, DEFAULT='linear')
%
% OUTPUT:
%  imOut        - Resampled image matrix
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Resamples an input image imIn oriented according to the homogeneous matrix matIn to and output
%              image imOut that has dimension dimOut and matrix matOut. This allows to resample images between
%              two spaces with different orientation and matrix sizes. It uses the Matlab interp3 function and
%              the interpolation method used by this method can be chosen. Note that the recommended use of this
%              function is resampling between space with similar resolution or upsampling. For downsampling, simple
%              interpolation does not delivery correct results respecting the point-spread-function and a combination of
%              xASL_im_PreSmooth and xASL_spm_reslice should be used instead.
%
% EXAMPLE:     imOut = xASL_im_ResampleIM(imIn,[1 0 0 10;0 1 0 2;0 0 1 3; 0 0 0 1],[31 0 0 0;0 2 0 2;0 0 2 3; 0 0 0 1],[50 50 50], 'linear')
%
% __________________________________
% Copyright 2015-2021 ExploreASL

% Input parameter admin
if nargin < 4
	error('Need to specify input image, dimensions and input and output orientation matrix');
end

if nargin < 5 || isempty('interpolationType')
    % Default, CAT12 also uses linear interpolation & single precision:
    interpolationType  = 'linear';
end

% Dimension of the input
dimIn = size(imIn);

% Define the coordinates of the axis of the input image
Xin = 1:(dimIn(2)-0);
Yin = 1:(dimIn(1)-0);
Zin = 1:(dimIn(3)-0);

% Coordinates for all voxels in the output image
[Xout,Yout,Zout]    = meshgrid(1:(dimOut(2)-0),1:(dimOut(1)-0),1:(dimOut(3)-0));

Xout = Xout(:);
Yout = Yout(:);
Zout = Zout(:);

% Calculate the coordinates in the input image
XYZresampled = (matIn\matOut) * [Xout';Yout';Zout';ones(size(Zout'))];

% Run the interpolation
imOut = interp3(Xin',Yin',Zin',double(imIn),XYZresampled(1,:),XYZresampled(2,:),XYZresampled(3,:),interpolationType);
imOut(isnan(imOut)) = 0; % remove nans
imOut = reshape(imOut,dimOut); % reshape back

end
