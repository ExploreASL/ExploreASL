function [imOut] = xASL_im_ResampleIM(imIn,matIn,matOut,dimOut,intMethod)
% xASL_im_ResampleIM Resamples an image using Matlab interp3 function
%
% IntMethod can be one of 'linear' 'cubic' 'or 'nearest'
%
% FORMAT: [imOut] = xASL_im_ResampleIM(imIn,matIn,matOut,dimOut,intMethod)
%
% INPUT:
%   imIn	    - ...
%   matIn	    - ...
%   matOut	    - ...
%   dimOut	    - ...
%   intMethod   - ...
%
% OUTPUT:
%  imOut        - ...
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Resamples an image using Matlab interp3 function.
%
% EXAMPLE:     n/a
%
% __________________________________
% Copyright 2015-2021 ExploreASL

if ~exist('intMethod','var')
    % default, CAT12 also uses linear interpolation & single precision:
    intMethod  = 'linear';
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
imOut = interp3(Xin',Yin',Zin',double(imIn),XYZresampled(1,:),XYZresampled(2,:),XYZresampled(3,:),intMethod);
imOut(isnan(imOut)) = 0; % remove nans
imOut = reshape(imOut,dimOut); % reshape back

return

