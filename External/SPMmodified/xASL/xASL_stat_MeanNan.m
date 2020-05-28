function y = xASL_stat_MeanNan(x,dim)
% Calculates the mean of the input while ignorning NaNs.
%
% FORMAT: y = xASL_stat_MeanNan(x[,dim])
%
% INPUT:
%   x 	   - input vector/matrix
%   dim    - dimension along which it operates (DEFAULT first non-singleton or 1)
% OUTPUT:
%   y      - calculated mean
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: It calculates the sum using the SUM functions and divides by the number of values but ignoring NaNs.
%
% EXAMPLE: y = xASL_stat_MeanNan([2 3 1; 0 -1 3])
%          y = xASL_stat_MeanNan([2 3 1; 0 -1 3],1)
%          y = xASL_stat_MeanNan([2 3 1; 0 -1 3],2)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright © 2015-2019 ExploreASL
%
% 2017-00-00 JP

x = double(x); % Single failed in large arrays according to CAT12

% Finds all the nans and non-nans
index_nans = isnan(x);
index_normals = (~index_nans);

% Sets nans to zero in the original data
x(index_nans) = 0;

% If working along the first dimension
if nargin < 2 || isempty(dim)
    % Number of the normal numbers
    count_normals = sum(index_normals);
    
    % Sum the numbers
    y = sum(x);
else
    % Count up non-NaNs.
    count_normals = sum(index_normals,dim);
    
    % Sum the numbers
    y = sum(x,dim);
end

% Set the lines with no non-nans to NaN to avoid division by 0
count_normals(count_normals==0) = NaN;

% Calculate the mean by dividing the sum with number of non-nans
y = y ./ count_normals;
end
