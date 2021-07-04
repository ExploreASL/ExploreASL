function y = xASL_stat_QuantileNan(x,quant,dim)
% Calculates the quantile of the input while ignorning NaNs.
%
% FORMAT: y = xASL_stat_QuantileNan(x[,quant, dim])
%
% INPUT:
%   x 	   - input vector/matrix
%   quant  - quantile between 0 and 1 (DEFAULT 0.5)
%   dim    - dimension along which it operates (DEFAULT first non-singleton or 1)
% OUTPUT:
%   y      - calculated quantile
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Calculates a quantile, but ignoring NaNs in the calculation
%
% EXAMPLE: y = xASL_stat_QuantileNan([2 3 1; 0 -1 3])
%          y = xASL_stat_QuantileNan([2 3 1; 0 -1 3],[],1)
%          y = xASL_stat_QuantileNan([2 3 1; 0 -1 3],0.5,2)
%          y = xASL_stat_QuantileNan([2 3 1; 0 -1 3],0.5)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright (c) 2015-2021 ExploreASL
%
% 2017-00-00 JP

%% Check input arguments
if nargin < 2 || isempty(quant)
	quant = 0.5;
end

% quant needs to be between 0 and 1 if outside this range, then throw and error
if (quant < 0) || (quant > 1)
	error('xASL_stat_QuantileNan: the value needs to be between 0 and 1');
end

%% Determine result
if isempty(x)
	if nargin < 3
		
		% The output size for [] is a special case when DIM is not given.
		if isequal(x,[])
			if isinteger(x)
				y = zeros(class(x));
			else
				y = nan(class(x));
			end
			return;
		end
		
		% Determine first nonsingleton dimension
		dim = find(size(x)~=1,1);
		
	end
	
	s = size(x);
	if dim <= length(s)
		s(dim) = 1;                  % Set size to 1 along dimension
	end
	if isinteger(x)
		y = zeros(s,class(x));
	else
		y = nan(s,class(x));
	end
	
elseif nargin < 3 && isvector(x)
	% If input is a vector, calculate single value of output.
	x = double(x); % Single failed in large arrays according to CAT12
	
	x = sort(x);
	nCompare = numel(x)-sum(isnan(x));
	half = getPosQuantile(quant,nCompare);
	y = x(half);
	if isnan(x(nCompare))        % Check last index for NaN
		y = nan(class(x));
	end
else
	x = double(x); % Single failed in large arrays according to CAT12
	if nargin < 3              % Determine first nonsingleton dimension
		dim = find(size(x)~=1,1);
	end
	
	s = size(x);
	
	if dim > length(s)           % If dimension is too high, just return input.
		y = x;
		return
	end
	
	% Permute the dimensions so that the dimension we work along is the first
	if dim > 1
		permuteDims = [dim:length(s) 1:(dim-1)];
		x = permute(x,permuteDims);
	else
		permuteDims = [];
	end
	
	snew = size(x);
	snewcols = prod(snew(2:end));
	% If more than 2D than reshape for 2D
	if length(snew) > 2
		x = reshape(x,snew(1), snewcols);
	end
	
	% Calculate the dimension after median in the new order
	sout = snew;
	sout(1) = 1;
	
	% Sort along the first dimension
	x = sort(x);
	
	% Indices of the nans
	nans = isnan(x);
	
	% Number of nonnans for each column
	nnonnans = snew(1)-sum(nans,1);
	
	% There are nans, we have to calculate the medians columnwise
	if sum(nans(:))
		y = nan(1,snewcols,class(x));
		% For each column
		for j = 1:snewcols
			% If there are non-nans
			if nnonnans(j)
				half = getPosQuantile(quant,nnonnans(j));
				y(:,j) = x(half,j);
			end
		end
	else % There are no nans
		nCompare = s(dim);           % Number of elements used to generate a median
		half = getPosQuantile(quant,nCompare);
		
		% If calculating along columns, use vectorized method with column
		% indexing.  Reshape at end to appropriate dimension.
		y = x(half,:);
		
		y(isnan(x(nCompare,:))) = NaN;   % Check last index for NaN
	end
	
	% Reshape back
	y = reshape(y,sout);
	
	% Permute the dimensions back
	if ~isempty(permuteDims)
		y = ipermute(y,permuteDims);
	end
	
end
end

% calculate the position quant from N
function n=getPosQuantile(quant,N)
n = round(N*quant);
if n<1
	n = 1;
end

if n>N
	n = N;
end

end