function new_mask = xASL_im_DilateErodeFull(mask,type,kernel)
% Runs dilation or erosion on a binary mask in full three dimensions
% It uses its own dilate_erode function and crops the image so that it
% contains only the mask

% Works only with odd sized kernels

% PARAMETERS
% mask = binary mask
% type == 'dilate' or 'erode'
% kernel - 3D vectors NxMxO containing the separable erosion kernel -
% can have different sizes.
new_mask = mask;

if (mod(size(kernel,1),2)+mod(size(kernel,2),2)+mod(size(kernel,3),2))<3
    error('Only odd kernel sizes accepted');
end

% If the input mask is empty then return an empty mask
if ~any(mask(:))
    return;
end

kernel_radius = floor(size(kernel)/2);
if length(kernel_radius) == 2
    kernel_radius(3) = 0;
end

% Sums around one of the dimension to discover the non-zero values 
vec1 = (1:size(mask,1))';
sum1 = squeeze(sum(mask,1)>0);

vec2 = (1:size(mask,2))';
vec3 = (1:size(mask,3))';
sum3 = sum(mask,3)>0;

% 3D version
[vecX,vecOrigX,vecBackX,vecY,vecOrigY,vecBackY,vecZ,vecOrigZ,vecBackZ,nonzero] = check_max3d((sum(sum3,2)>0).*vec1,(sum(sum1,2)>0).*vec2,(sum(sum1,1)>0).*vec3',size(mask,1),size(mask,2),size(mask,3),kernel_radius);
if nonzero
    % Do the padding there and back - only once for all dimensions
    % Run the separable dilation - specify the vector - always Nx1 and
    % specify the dimension - all three will have separate code because the
    % for will run across different voxels
    
    % Cropping or padding of the mask
    mask_cropped = mask(vecX,vecY,vecZ);
    
    % Run the dilation/erosion through mex files - separably in all three
    % dimension
    mask_cropped = xASL_mex_dilate_erode_3D(double(mask_cropped),double(kernel),type);
    
    % Bring the processed image to the original space 
    new_mask(vecOrigX,vecOrigY,vecOrigZ) = mask_cropped(vecBackX,vecBackY,vecBackZ);
    
end
return;

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

return

