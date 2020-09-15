%cat_vol_median3c Median filter for label maps.
%  Median Filter for 3D single image D. Bi is used to mask voxels for the 
%  filter process, whereas Bn is used to mask voxels that are used as 
%  neighbors in the filter process. 
% 
%  This is a special subversion to filter label maps!
% 
%  M = cat_vol_median3c(D[,Bi,Bn])
% 
%  D  (single)  .. 3D matrix for filter process 
%  Bi (logical) .. 3D matrix that mark voxel that should be filtered
%  Bn (logical) .. 3D matrix that mark voxel that are used as neighbors 
% 
%  See also cat_vol_median3, compile.
%  ________________________________________________________________________
%  Robert Dahnke, 2012/01 
%  $Id: cat_vol_median3c.m 1523 2019-11-21 23:12:24Z gaser $
