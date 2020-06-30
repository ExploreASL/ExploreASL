function [mu,su,nu] = kmeans3D(y,k)
% K-means clustering
%_______________________________________________________________________
% FORMAT [mu,su,nu] = kmeans3D(y,k)
% 
%  y  .. data 
%  k  .. Number of components
%
%  mu .. vector of class means 
%  su .. vector of class std 
%  nu .. vector of class number of voxel
%
% modified version of
% spm_kmeans1.m 1143 2008-02-07 19:33:33Z spm $
%_______________________________________________________________________
% Christian Gaser & Robert Dahnke
% $Id: kmeans3D.m 1368 2018-09-28 10:04:52Z dahnke $

if nargin<1, help kmeans3D; return; end
if nargin<2, k=1; end

k = max(1,k); 

y = y(:)';
y(isnan(y))=[]; % remove NaNs 
if numel(y)<=0
  mu = nan(1,k);
  su = mu; 
  nu = mu; 
  return; 
end

N = length(y);

% Spread seeds evenly according to CDF
x = sort(y);
seeds=[1,2*ones(1,k-1)]*N/(2*k);
seeds=ceil(cumsum(seeds));

last_i = N; %(ones(1,N);
mu = x(seeds);
su = zeros(size(mu));
nu = ones(size(mu));

d = zeros(k,length(y));
for loops = 1:1000,  
  for j=1:k,
    d(j,:) = (y-mu(j)).^2;
  end
  [tmp,i] = min(d); clear tmp %#ok<ASGLU>
  if sum(i - last_i)==0 || isempty(y(i==j))
    % If assignment is unchanged
    break;
  else
   % Recompute centres
   for j=1:k,
     mu(j) = mean(y(i==j));
     su(j) = std(y(i==j));
     nu(j) = sum(i==j)./numel(y(:));
   end
   last_i=i;
  end
end  
