%cat_vol_ROIval Region-wise statistic.
%  Estimation of mean, standard deviation, minimum, maximum, sum, number
%  of Yv in a ROI described by an label map Ya.
%
%  [mn,std,min,max,sum,num] = cat_vol_ROIval(Ya,Yv)
%
%  mn   (single)    .. mean value
%  std  (single)    .. standard deviation 
%  min  (single)    .. minimum of Yv each ROI
%  max  (single)    .. maximum of each ROI
%  sum  (single)    .. sum of all values of each ROI
%  num  (single)    .. number of voxel of each ROI
%  Ya   (3D uint8)  .. label volume
%  Yv   (3D single) .. data volume
% 
%  See also compile.
%  ________________________________________________________________________
%  $Id: cat_vol_ROIval.m 1531 2019-12-05 16:07:30Z gaser $
