%cat_vol_eidist Voxel-based eiknoal distance calculation.
%  This function estimates the Euclidean distance D to the closest boundary 
%  voxel I given by the 0.5 isolevel in B that should contain values 
%  between 0 and 1 to define the boundary using PVE.  To align the closest 
%  voxel a modified Eikonal distances is estimated on the field L.
%
%  For a correct side alignment a harder option "csf" is used that avoids
%  the diffusion to voxels with greater values of L (L[i]>L[ni] for
%  diffusion).
%
%  Voxels with NaN and -INF are ignored and produce NaN in D, whereas 
%  in I the default index is used, because I has to be a integer matrix 
%  where NaN is not available. With setnan=0, the result of NAN voxel in 
%  D is changed to INF.
% 
%  [D,I,Dw] = cat_vol_eidist(B,L,[vx_vol,euclid,csf,setnan,verb])
% 
%  D      .. Euclidean distance map to the nearest Boundary point in the
%             Eikonal field (3d-single-matrix)
%            (air distance)
%  D      .. Euclidean distance map in the Eikonal field (3d-single-matrix)
%            (way length) 
%  I      .. index map      (3d-uint32-matrix)
%  B      .. boundary map   (3d-single-matrix)
%  L      .. speed map      (3d-single-matrix)
%  vx_vol .. voxel size     (1x3-double-matrix): default=[1 1 1]
%            ... not tested yet
%  csf    .. option         (1x1-double-value): 0-no; 1-yes; default=1
%  euclid .. option         (1x1-double-value): 0-no; 1-yes; default=1
%            output euclidean or speed map 
%  setnan .. option         (1x1-double-value): 0-no; 1-yes; default=1
%  verb   .. option         (1x1-double-value): 0-no, 1-yes; default=0
% 
%
%  Examples: 
%  Definitions of a object matrix A and of a speedmap F:
%
%  1) 
%    A = zeros(50,50,50,'single'); 
%    A(20:30,5:15,20:30) = 1; A(20:30,35:45,20:30) = 1; 
%    A = smooth3(A); A(1:5,1:25,:) = nan; 
%    A(1:5,26:50,:) = -inf; A(45:50,26:50,:) = inf;
%    F = ones(size(A),'single'); F(10:40,20,:) = 0.5; F(40,10:40,:) = 0;
%
%  2) 
%    A = zeros(10,20,10,'single'); 
%    A(:,1:5,:) = 1; A(:,15:end,:) = nan; 
%    F = ones(size(A),'single');
%
%  3)
%    A = zeros(10,20,10,'single'); 
%    A(:,1:5,:) = 1; A(:,6,:) = 0.2; A(:,15:end,:) = nan; 
%    F = ones(size(A),'single');
%
%  4)
%    A = zeros(10,20,10,'single'); A(:,1:5,:) = 1; A(:,15:end,:) = nan; 
%    F = ones(size(A),'single');
% 
%  Process and show data:
%    [D,I] = cat_vol_eidist(A,F,[1 1 1],1,1,0,1);
%    ds('x2','',1,A,D,D,I,round(size(A,3)/3));
%
%  See also cat_vbdist, compile, ds.
%  ________________________________________________________________________
%  Robert Dahnke 2010/01
%  $Id: cat_vol_eidist.m 1523 2019-11-21 23:12:24Z gaser $
