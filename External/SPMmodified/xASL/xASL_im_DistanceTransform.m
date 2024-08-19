function [dist, x, y, z] = xASL_im_DistanceTransform(im)
% xASL_im_DistanceTransform Calculates the distance transform in a binary image
% Uses Borgefors Chamfers computation of Euclidean distance in 3D using a 5x5x5 window
%
% FORMAT: [dist, x, y, z] = xASL_im_DistanceTransform(im)
%
% INPUT:
%   im	    - input image
%
% OUTPUT:
%  dist     - distance map, distance in voxels to the closest point on the mask
%  x,y,z    - X, Y, and Z coordinate of the closest voxel on the mask. Note that these are relative coordinates describing the relative coordinates from the current voxel
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Calculates the distance transform in a binary image
% Uses Borgefors Chamfers computation of Euclidean distance in 3D using a
% 5x5x5 window.
%
% EXAMPLE:     n/a
%
% Implemented by Jan Petr
% According to: Stina Svensson and Gunilla Borgefors, Computer Vision and Image Understanding 88, 24–53 (2002)
% doi:10.1006/cviu.2002.0976
%
% __________________________________
% Copyright 2015-2024 ExploreASL
 
imBin = double( im > 0);

% Note the the MEX file CHAMFERS3D have Y and X reversed with respect to the Matlab notation, so having y,x,z is correct here.
[dist, y, x, z] = xASL_mex_chamfers3D(imBin);

return;
