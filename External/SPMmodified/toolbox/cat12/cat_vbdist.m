%cat_vbdist Voxel-based euclidean distance calculation.
%  Calculates the euclidean distance without partial volume effect to an 
%  object in P with a boundary of 0.5.
% 
%  [D,I,L] = cat_vbdist(P[,R])
%  
%  D (single)  .. distance image
%  I (uint32)  .. index of nearest point
%  L (uint8)   .. label map
%  P (single)  .. input image 
%  R (logical) .. range for distance calculation (to speedup processing)
%
%  See also bwdist, cat_vol_eidist, compile.
%  ________________________________________________________________________
%  Robert Dahnke 2010/01
%  $Id: cat_vbdist.m 1523 2019-11-21 23:12:24Z gaser $
