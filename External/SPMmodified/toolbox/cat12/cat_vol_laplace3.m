%cat_vol_laplace3 Volumetric Laplace filter with Dirichlet boundary. 
%  Filter SEG within the intensity range of low and high until the changes
%  are below TH. 
% 
%  L = cat_vol_laplace3(SEG,low,high,TH)
%
%  SEG  .. 3D single input matrix
%  low  .. low boundary threshold
%  high .. high boundary threshold
%  TH   .. threshold to control the number of iterations
%          maximum change of an element after iteration
%
%  See also cat_vol_laplace3R, compile.
%  ________________________________________________________________________
%  Robert Dahnke 2009/01
%  $Id: cat_vol_laplace3.m 1561 2020-02-04 15:49:34Z gaser $
