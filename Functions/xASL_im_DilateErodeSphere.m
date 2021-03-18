function el = xASL_im_DilateErodeSphere(R)
%xASL_im_DilateErodeSphere Return a 3D structuring element (binary) sphere of a given diameter
%
% FORMAT:       el = xASL_im_DilateErodeSphere(R)
% 
% INPUT:        R - diameter of the sphere (OPTIONAL, DEFAULT = 1)
%
% OUTPUT:       el - the structural element
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Creates a 3D structuring element (binary) sphere with the given diameter (R) and size 2*R+1
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:  el = xASL_im_DilateErodeSphere(2)
% __________________________________
% Copyright 2015-2020 ExploreASL

%% Admin
if nargin < 1 || isempty(R)
	R = 1;
end

% Size of the element
N=2*R+1;

% Field for distance calculation
[a,b,c] = meshgrid(1:N,1:N,1:N);

% Calculate the distance from the center
C=R+1;
d=sqrt((a-C).^2 + (b-C).^2 + (c-C).^2);

el = d<=R;

end
