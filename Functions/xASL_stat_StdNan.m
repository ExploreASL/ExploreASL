function y = xASL_stat_StdNan(varargin)
% Calculates standard deviation of values in X while ignoring NaNs.
%
% FORMAT: y = xASL_stat_StdNan(x[,w,dim])
%
% INPUT:
%   x 	   - input vector/matrix
%   w      - normalization/weighting (DEFAULT 0)
%            0 - normalizes by N-1 (N is the length of the dimension along which it operates)
%            1 - normalizes by N
%            W - computes the standard deviation using the weight vector W which has the length N
%   dim    - dimension along which it operates (DEFAULT first non-singleton or 1)
% OUTPUT:
%   y      - calculated standard deviation
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: It behaves in a similar way as VAR - it directly passes all arguments to xASL_stat_VarNan.
%
% EXAMPLE: y = xASL_stat_StdNan([2 3 1; 0 -1 3])
%          y = xASL_stat_StdNan([2 3 1; 0 -1 3],1)
%          y = xASL_stat_StdNan([2 3 1; 0 -1 3],0,2)
%          y = xASL_stat_StdNan([2 3 1; 0 -1 3],[1 1 3],2)
%          y = xASL_stat_StdNan([2 3 1; 0 -1 3],[],2)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright © 2015-2019 ExploreASL
%
% 2017-00-00 JP

y = sqrt(xASL_stat_VarNan(varargin{:}));
end
