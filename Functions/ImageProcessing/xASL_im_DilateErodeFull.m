function imOut = xASL_im_DilateErodeFull(imIn, type, kernel)
%xASL_im_DilateErodeFull Dilation or erosion of a binary imIn using a 3D kernel
%
% FORMAT: imOut = xASL_im_DilateErodeFull(imIn, type, kernel)
% 
% INPUT:  imIn   - binary input image (REQUIRED)
%         type   - a string indicating the operation type 'dilate' or 'erode' (REQUIRED)
%         kernel - 3D vector NxMxO containing 3D kernel
%
% OUTPUT: imOut - binary image after erosion/dilation
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Runs dilation or erosion on a binary imIn in full three dimensions.
%               It uses its own dilate_erode function and crops the image so that it
%               contains only the mask. The size of all three dimensions of the kernel needs to be an odd number.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: imOut = xASL_im_DilateErodeFull(imIn, 'erode', [1 1 1;0 1 0; 1 1 1])
% __________________________________
% Copyright 2015-2021 ExploreASL

%% Admin
if nargin < 1
	error('Input image needed');
end

if nargin < 2 || (~strcmpi(type,'dilate') && (~strcmpi(type,'erode')))
	error('Operation type has to be defined as dilate or erode');
end

if nargin < 3 || (mod(size(kernel,1),2)+mod(size(kernel,2),2)+mod(size(kernel,3),2))<3
    error('Only odd kernel sizes accepted');
end

imOut = imIn;

% If the input mask is empty then return an empty mask
if ~any(imIn(:))
    return;
end

%% Prepares the cropping
kernel_radius = floor(size(kernel)/2);
if length(kernel_radius) == 2
    kernel_radius(3) = 0;
end

% Sums around one of the dimension to discover the non-zero values 
vec1 = (1:size(imIn,1))';
sum1 = squeeze(sum(imIn,1)>0);

vec2 = (1:size(imIn,2))';
vec3 = (1:size(imIn,3))';
sum3 = sum(imIn,3)>0;

% 3D version
[vecX,vecOrigX,vecBackX,vecY,vecOrigY,vecBackY,vecZ,vecOrigZ,vecBackZ,nonzero] = check_max3d((sum(sum3,2)>0).*vec1,(sum(sum1,2)>0).*vec2,(sum(sum1,1)>0).*vec3',size(imIn,1),size(imIn,2),size(imIn,3),kernel_radius);

if nonzero
    % Cropping or padding of the mask
    mask_cropped = imIn(vecX,vecY,vecZ);
    
    % Run the dilation/erosion through mex files - separably in all three
    % dimension
    mask_cropped = xASL_mex_dilate_erode_3D(double(mask_cropped),double(kernel),type);
    
    % Bring the processed image to the original space 
    imOut(vecOrigX,vecOrigY,vecOrigZ) = mask_cropped(vecBackX,vecBackY,vecBackZ);
end

end

% Calculate for the given mask the cropped or enlarged image
% The advantage is that we enlarge the image so that we run the kernel in
% the middle and do not have to check for the boundary condition since the
% sides are padded
function [vecX,vecOrigX,vecBackX,vecY,vecOrigY,vecBackY,vecZ,vecOrigZ,vecBackZ,nonzero] = check_max3d(indX,indY,indZ,NX,NY,NZ,kernel_radius)
indX = indX(indX>0);
indY = indY(indY>0);
indZ = indZ(indZ>0);

if isempty(indX)
    minX = 0;
    maxX = 0;
else
    minX = min(indX);
    maxX = max(indX);
end

if isempty(indY)
    minY = 0;
    maxY = 0;
else
    minY = min(indY);
    maxY = max(indY);
end

if isempty(indZ)
    minZ = 0;
    maxZ = 0;
else
    minZ = min(indZ);
    maxZ = max(indZ);
end

if (maxX>0) && (maxY>0) && (maxZ > 0)
    nonzero = 1;
    minX = minX-2*kernel_radius(1);
    maxX = maxX+2*kernel_radius(1);
    vecX = [ones(1,-minX+1),max(1,minX):min(NX,maxX),NX*ones(1,maxX-NX)]; % Orig->padded working copy
    vecOrigX = max(1,minX):min(NX,maxX); % To these indexes is transfered the result image
    vecBackX = vecOrigX-minX+1; % These indexes from the result are copied to the indexes of the orginal above
        
    minY = minY-2*kernel_radius(2);
    maxY = maxY+2*kernel_radius(2);
    vecY = [ones(1,-minY+1),max(1,minY):min(NY,maxY),NY*ones(1,maxY-NY)];
    vecOrigY = max(1,minY):min(NY,maxY); % To these indexes is transfered the result image
    vecBackY = vecOrigY-minY+1; % These indexes from the result are copied to the indexes of the orginal above
    
    minZ = minZ-2*kernel_radius(3);
    maxZ = maxZ+2*kernel_radius(3);
    vecZ = [ones(1,-minZ+1),max(1,minZ):min(NZ,maxZ),NZ*ones(1,maxZ-NZ)];
    vecOrigZ = max(1,minZ):min(NZ,maxZ); % To these indexes is transfered the result image
    vecBackZ = vecOrigZ-minZ+1; % These indexes from the result are copied to the indexes of the orginal above
else
    nonzero = 0;
end

end

