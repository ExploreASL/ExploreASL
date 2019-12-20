function el = xASL_im_DilateErodeSphere(R)
% 3D structuring element (binary) sphere

% Size of the element
N=2*R+1;

% Field for distance calculation
[a,b,c] = meshgrid(1:N,1:N,1:N);

% Calculate the distance from the center
C=R+1;
d=sqrt((a-C).^2 + (b-C).^2 + (c-C).^2);

el = d<=R;

return;