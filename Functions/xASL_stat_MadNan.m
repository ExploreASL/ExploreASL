function y = xASL_stat_MadNan(x,flag,dim)
% Calculates the Median/mean absolute deviation. Use functions that ignore NANs.
%
% FORMAT: y = xASL_stat_MadNan(x[,flag, dim])
%
% INPUT:
%   x 	   - input vector/matrix
%   flag   - 0 Mean absolute deviation, 1 Median absolute deviation (DEFAULT 0)
%   dim    - dimension along which it operates (DEFAULT first non-singleton or 1)
% OUTPUT:
%   y      - calculated MAD value
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Calculates a Median/Mean Absolute deviation, but ignoring NaNs in the calculation.
%   xASL_stat_MadNan(X) or xASL_stat_MadNan(X,0) computes xASL_stat_MeanNan(ABS(X-xASL_stat_MeanNan(X))
%   xASL_stat_MadNan(X,1) computes xASL_stat_MedianNan(ABS(X-xASL_st_MedianNan(X)).
%
%
% EXAMPLE: y = xASL_stat_MadNan([2 3 1; 0 -1 3],0)
%          y = xASL_stat_MadNan([2 3 1; 0 -1 3],1,1)
%          y = xASL_stat_MadNan([2 3 1; 0 -1 3],[],2)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright © 2015-2019 ExploreASL
%
% 2017-00-00 JP

% Set the flag if missing
if nargin < 2 || isempty(flag)
    flag = 0;
end

% For empty input return Nan
if isequal(x,[])
    y = NaN;
    return;
end

x = double(x); % Single failed in large arrays according to CAT12

% Set the dimension if missing
if nargin < 3 || isempty(dim)
	% Figure out which dimension nanmean will work along.
	dim = find(size(x) ~= 1, 1);
	if isempty(dim)
		dim = 1;
	end
end

% Need to tile the output of nanmean to center X.
rpsize = ones(1,max(ndims(x),dim));
rpsize(dim) = size(x,dim);

if flag == 0
    % Mean absolute deviation
    y = xASL_stat_MeanNan( abs(x - repmat(xASL_stat_MeanNan(x,dim), rpsize)) ,dim);
else
    % Median absolute deviation
    y = xASL_stat_MedianNan( abs(x - repmat(xASL_stat_MedianNan(x,dim), rpsize)) ,dim);
end
    
end
