function bIsnear = xASL_stat_IsEqualTol(x, y, tol)
% Verifies if two numerical arrays are nearly equal
%
% FORMAT: bIsnear = xASL_stat_IsEqualTol(x, y[, tol])
%
% INPUT:
%   x 	    - the first input scalar/vector/matrix (REQUIRED)
%   y 	    - the second input scalar/vector/matrix (REQUIRED)
%   tol     - tolerance scalar/vector/matrix (OPTIONAL, DEFAULT 1e-8)
% OUTPUT:
%   bIsnear - output, true if values are nearly equal
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Checks if two input numerical arrays X and Y are equal or nearly equal with a given tolerance TOL. The output is a logical array
%              that is TRUE in all positions where X and Y are nearly equal. X, Y, tol can be all scalars or all matrices with the same size. Any combination
%              of scalars and matrices can be given on the input given that the matrices are equal in size.
%
% EXAMPLE: 
%    bIsnear = xASL_stat_IsEqualTol(1, 2, 3)
%    bIsnear = xASL_stat_IsEqualTol(1, [2,3,4], 0.1)
%    bIsnear = xASL_stat_IsEqualTol([1,2;3 4], [1,2;3,4], 0.1)
%    bIsnear = xASL_stat_IsEqualTol(1, 2)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright 2015-2024 ExploreASL


% Admin
if nargin < 2 || isempty(x) || isempty(y)
    error('Two non-empty inputs are required')
end

% Default tolerance
if nargin < 3 || isempty(tol)
   tol = 1e-8;
end

% Check that all inputs are numerical arrays
if ~isnumeric(x) || ~isnumeric(y) || ~isnumeric(tol)
   error('All inputs must be numerical arrays')
end

% Check if X and Y inputs are the same size or scalars
if ~isequal(size(x), size(y)) && numel(x) > 1 && numel(y) > 1
   error('Inputs X and Y must be the same size or either of them must be a scalar');
end

% Check if tol is a scalar or has the same size as X and Y
if numel(tol) > 1 && ( (numel(x) > 1 && ~isequal(size(tol), size(x))) || (numel(y) > 1 && ~isequal(size(tol), size(y))) )
	error('Input TOL must have the same size as X or Y or be a scalar');
end

% Check if the difference is smaller than tolerance
bIsnear = abs(x-y) <= abs(tol);

end
