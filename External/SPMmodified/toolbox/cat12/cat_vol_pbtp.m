%cat_vol_pbtp Projection-Based Thickness.
%  Main thickness projection function (<a href="matlab:web('http://dbm.neuro.uni-jena.de/pdf-files/Dahnke-NI12.pdf','-browser');">Dahnke et al., 2011)</a>).
%
%  [GMT,RPM] = cat_vol_pbtp(SEG,WMD,CSFD[,opt])
%
%  GMT   (3D  single) .. thickness image
%  RPM   (3D  single) .. radial position map
%  SEG   (3D  single) .. segment image with low and high boundary bd
%                        (default 1=CSF, 2=GM, 3=WM)
%  WMD   (3D  single) .. CSF distance map
%  CSFD  (3D  single) .. CSF distance map
%
%  opt                .. MATLAB structure
%   .bd  (1x2 single) .. [low,high] boundary values (default 1.5 and 2.5)
%   .CSFD             .. calculate CSFD
%   .PVE              .. use PVE information (0=none,1=fast,2=exact)
%
%  See also cat_vbdist, cat_vol_eidist, compile.
%  ________________________________________________________________________
%  Robert Dahnke 2009/11, 2011/12
%  $Id: cat_vol_pbtp.m 1524 2019-11-28 17:00:31Z gaser $
