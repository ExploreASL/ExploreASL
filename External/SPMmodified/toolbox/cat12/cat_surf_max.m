function [N,Z,M,A,XYZ] = cat_surf_max(X,L,G)
% Sizes, local maxima and locations of excursion sets on a surface mesh
% FORMAT [N,Z,M,A,XYZ] = cat_surf_max(X,L,G)
% X        - a [nx1] array of stat values
% L        - a [nx1] array of locations {in vertices}
% G        - a patch structure
%
% N        - a [px1] size of connected components {in vertices}
% Z        - stat values of maxima
% M        - location of maxima {in vertices}
% A        - region number
% XYZ      - cell array of vertices locations
%__________________________________________________________________________
%
% See also: spm_max.m, spm_mesh_clusters.m and spm_mesh_get_lm.m
%__________________________________________________________________________
% Copyright (C) 2012-2016 Wellcome Trust Centre for Neuroimaging
%
% modified version of spm_mesh_max because for some rare data only one global 
% maximum was returned
% Guillaume Flandin
% $Id: cat_surf_max.m 1281 2018-02-18 22:05:09Z gaser $
%
% $Id: cat_surf_max.m 1281 2018-02-18 22:05:09Z gaser $

%-Get connected components
%--------------------------------------------------------------------------
LL     = NaN(size(G.vertices,1),1);
LL(L(1,:)) = X;
[C, N0] = spm_mesh_clusters(G,LL);

%-Get local maxima
%--------------------------------------------------------------------------
M0 = spm_mesh_get_lm(G,LL);

Z = LL(M0);
A = C(M0);
M = [M0;ones(2,size(M0,2))];
N = N0(A);

% In rare cases where local maxima are not correctly found the number of clusters 
% differs from cluster indices in A. Then, we have to add non-unique entries 
% using function find_connected_component
if numel(unique(N0)) ~= max(A)
  A2 = spm_mesh_adjacency(G);
  A2 = A2 + speye(size(A2));
  M2 = find_connected_component(A2,LL);
  
  % find uniqe entries in both methods
  M = [M0 M2];
  [tmp, ind] = unique(M);
  M0 = M(sort(ind));
  
  % and create new variables using corrected M0
  Z = LL(M0);
  A = C(M0);
  M = [M0;ones(2,size(M0,2))];
  N = N0(A);
end

if nargout > 4
  XYZ = cell(1,max(A));
  for i=1:numel(XYZ)
    XYZ{i} = find(C==i)';
    XYZ{i} = [XYZ{i};ones(2,size(XYZ{i},2))];
  end
end

%==========================================================================
function C = find_connected_component(A, T);
% find connected components 
% FORMAT C = find_connected_component(A,T)
% A        - a [nxn[ (reduced) adjacency matrix
% T        - a [nx1] data vector (using NaNs or logicals), n = #vertices
%
% C        - a [nx1] vector of cluster indices
%
% modified version from spm_mesh_clusters.m 5065 2012-11-16 20:00:21Z guillaume
%

%-Input parameters
%--------------------------------------------------------------------------
T0 = T;
if ~islogical(T)
  T   = ~isnan(T);
end
  
A1 = A;
A1(~T,:) = [];
A1(:,~T) = [];

%-And perform Dulmage-Mendelsohn decomposition to find connected components
%--------------------------------------------------------------------------
[p,q,r] = dmperm(A1);
N       = diff(r);
CC      = zeros(size(A1,1),1);
for i = 1:length(r)-1
  CC(p(r(i):r(i+1)-1)) = i;
end
C       = NaN(numel(T),1);
C(T)    = CC;

mx_array  = zeros(1,max(C));
ind_array = zeros(1,max(C));

for i = 1:max(C)
  N = find(C == i);
  T1 = zeros(size(T0));
  T1(N) = T0(N);
  [mx_array(i), ind_array(i)] = max(T1);
end

%-Sort connected component labels according to their max T value
%--------------------------------------------------------------------------
[tmp,ni]  = sort(mx_array, 'descend');
C = ind_array(ni);
