%cat_vol_simgrow Volumetric region-growing.
%  
%  [SLAB,DIST] = cat_vol_simgrow(ALAB,SEG[,d,dims,dd]);
%  
%  SLAB (3D  single) .. output label map  
%  DIST (3D  single) .. distance map from region-growing
%  ALAB (3D  single) .. input label map
%  SEG  (3D  single) .. input tissue map
%  d    (1x1 double) .. growing treshhold paramter (max local gradient)
%                       in SEG
%  dims (1x3 double) .. voxel dimensions (default [1,1,1])
%  dd   (1x2 double) .. general growing limits in SEG 
%
%  See also cat_vol_downcut.
%  ________________________________________________________________________
%  Robert Dahnke 2010/01
%  $Id: cat_vol_simgrow.m 1524 2019-11-28 17:00:31Z gaser $
