function y = xASL_stat_SumNan(x,dim)
% Calculates the sum of the input while ignorning NaNs.
%
% FORMAT: y = xASL_stat_SumNan(x[,dim])
%
% INPUT:
%   x 	   - input vector/matrix
%   dim    - dimension along which it operates (DEFAULT first non-singleton or 1)
% OUTPUT:
%   y      - calculated sum
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: It uses the function SUM, but it sets all the NaNs to zero before calling it.
%
% EXAMPLE: y = xASL_stat_SumNan([2 3 1; 0 -1 3])
%          y = xASL_stat_SumNan([2 3 1; 0 -1 3],1)
%          y = xASL_stat_SumNan([2 3 1; 0 -1 3],2)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright © 2015-2019 ExploreASL
%
% 2017-00-00 JP

x = double(x); % Single failed in large arrays according to CAT12

% Sets nans to zero
x(isnan(x)) = 0;

% if the dimension is not given
if nargin < 2  
    y = sum(x);
% if the dimension is given
else           
    y = sum(x,dim);
end

end