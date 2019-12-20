function rotated = xASL_im_rotate( im, angle )
% Simple rotation of the first two dimension of a ND image by 0,90,180,270 degrees
% INPUT:
% 	im - an ND image
%	angle - angle of rotation in the first two dimensions, can be either 0,90,180,270 or other multiple of 90
% OUTPUT:
%	rotated - rotated image of the same dimensions
% JP, xASL, 2018

% In case the angle is higher than 360 or smaller than 0
angle = mod(angle,360);

% Skip empty images
if (isempty(im)) 
    rotated = im;
    return;
end

% If the number of dimension is higher than 3, then reshape to 3D
orig_size = size(im);
new_size = orig_size;
is4D = 0;
if ndims(im) > 3
	im = reshape(im,orig_size(1),orig_size(2),[]);
	is4D = 1;
end

switch (angle)
	% Do nothing for 0 degrees
    case 0
        rotated = im;
    case 90
		% Calculate the size of the new image
		new_size(1) = orig_size(2);
		new_size(2) = orig_size(1);

		% Perform the 90 degree rotation
        rotated = zeros([size(im,2),size(im,1),size(im,3)]);
		for k = 1:size(im,3)
			rotated(:,:,k) = im(:,:,k).';
			rotated(:,:,k) = rotated(end:-1:1,:,k);
		end
        
		% Reshape the ND file
		if is4D
			rotated = reshape(rotated,new_size);
		end
    case 180
        rotated = reshape(im(end:-1:1,end:-1:1,:),size(im));
    case 270
		% Calculate the size of the new image
		new_size(1) = orig_size(2);
		new_size(2) = orig_size(1);

       % Perform the 270 degree rotation
        rotated = zeros([size(im,2),size(im,1),size(im,3)]);
		for k = 1:size(im,3)
            rotated(:,:,k) = im(:,:,k).';
            rotated(:,:,k) = rotated(:,end:-1:1,k);
		end
        
		% Reshape the ND file
		if is4D
			rotated = reshape(rotated,new_size);
		end
    otherwise
        error('xASL_im_rotate: Invalid input, only does rotation by 0, 90, 180 or 270 degrees.');
end

return;

