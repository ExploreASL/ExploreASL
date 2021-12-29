function y = xASL_stat_VarNan(x,w,dim)
% Calculates variance of values in X while ignoring NaNs.
%
% FORMAT: y = xASL_stat_VarNan(x[,w,dim])
%
% INPUT:
%   x 	   - input vector/matrix
%   w      - normalization/weighting (DEFAULT 0)
%            0 - normalizes by N-1 (N is the length of the dimension along which it operates)
%            1 - normalizes by N
%            W - computes the variance using the weight vector W which has the length N
%   dim    - dimension along which it operates (DEFAULT first non-singleton or 1)
% OUTPUT:
%   y      - calculated variance
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: It behaves in a similar way as VAR.
%
% EXAMPLE: y = xASL_stat_VarNan([2 3 1; 0 -1 3])
%          y = xASL_stat_VarNan([2 3 1; 0 -1 3],1)
%          y = xASL_stat_VarNan([2 3 1; 0 -1 3],0,2)
%          y = xASL_stat_VarNan([2 3 1; 0 -1 3],[1 1 3],2)
%          y = xASL_stat_VarNan([2 3 1; 0 -1 3],[],2)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright © 2015-2019 ExploreASL
%
% 2017-00-00 JP


% Admin
if nargin < 2 || isempty(w)
	w = 0; 
end

if (nargin < 3 || isempty(dim)) && w>1 %% UGLY QUICK FIX
                                        % allows to input dimension along
                                        % which we operate, as second
                                        % argument, similar to other
                                        % operations such as mean, sum etc
                                        % Allowing to iterate over
                                        % operations
    dim = w;
    w = 0;    
elseif (nargin < 3 || isempty(dim))
	% The output size for [] is a special case when DIM is not given.
	if isequal(x,[])
		y = NaN(class(x));
		return;
	end
	
	% Figure out which dimension sum will work along - takes first non-singleton
	dim = find(size(x) ~= 1, 1);
	if isempty(dim)
		dim = 1;
	end
end

x = double(x); % Single failed in large arrays according to CAT12

% Unweighted variance
if isequal(w,0) || isequal(w,1)
    n = sum(~isnan(x),dim);
    
    if w == 0 
        % The unbiased estimator: divide by (n-1).  Can't do this
        % when n == 0 or 1, so for n == 1 the result will be NaN
        denom = max(n - 1,1);
    else
        % The biased estimator: divide by n.
        denom = n; % n==0 => return NaNs, n==1 => return zeros
    end
    denom(n==0) = NaN;

    xmean = xASL_stat_MeanNan(x, dim);
    x = bsxfun(@minus, x, xmean);

    y = xASL_stat_SumNan(abs(x).^2, dim) ./ denom; % abs guarantees a real result

    % Weighted variance
elseif isvector(w) && all(w >= 0)
    siz = size(x);
    if numel(w) ~= siz(dim)
        if isscalar(w)
            error(message('xASL_stat_VarNan:Invalid weights'));
        else
            error(message('xASL_stat_VarNan:Invalid size weights'));
        end
    end

    % Normalize W, and embed it in the right number of dims.  Then
    % replicate it out along the non-working dims to match X's size.
    wresize = ones(1,max(ndims(x),dim)); wresize(dim) = siz(dim);
    w = reshape(w, wresize);
    wrepmat = siz;wrepmat(dim) = 1;
    w = repmat(w,wrepmat);
    
    n = xASL_stat_SumNan(~isnan(x).*w,dim);
    
    x0 = xASL_stat_SumNan(x.*w,dim)./n;
    x = bsxfun(@minus, x, x0);
    y = xASL_stat_SumNan(bsxfun(@times, w, abs(x).^2), dim)./n; % abs guarantees a real result

else
    error(message('xASL_stat_VarNan:Invalid weights'));
end
end

