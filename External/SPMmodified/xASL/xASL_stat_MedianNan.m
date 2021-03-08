function y = xASL_stat_MedianNan(x,dim)
% Calculates the median of the input while ignorning NaNs.
%
% FORMAT: y = xASL_stat_MedianNan(x[,dim])
%
% INPUT:
%   x 	   - input vector/matrix
%   dim    - dimension along which it operates (DEFAULT first non-singleton or 1)
% OUTPUT:
%   y      - calculated median
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: It calculates the MEDIAN along the given dimension, but it sets all the NaNs to zero before calling it.
%
% EXAMPLE: y = xASL_stat_MedianNan([2 3 1; 0 -1 3])
%          y = xASL_stat_MedianNan([2 3 1; 0 -1 3],1)
%          y = xASL_stat_MedianNan([2 3 1; 0 -1 3],2)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright © 2015-2019 ExploreASL
%
% 2017-00-00 JP

if isempty(x)
  if nargin < 2 || isempty(dim)

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

elseif nargin < 2 && isvector(x)
  % If input is a vector, calculate single value of output.
  x = double(x); % Single failed in large arrays according to CAT12

  % Sort the vector
  x = sort(x);
  nCompare = numel(x)-sum(isnan(x));
  half = floor(nCompare/2);
  y = x(half+1);
  % If at least a single non-NaN element
  if nCompare > 0
	  if 2*half == nCompare        % Average if even number of elements
		  y = (x(half)+y)/2;
	  end
	  if isnan(x(nCompare))        % Check last index for NaN
		  y = nan(class(x));
	  end
  else
	  % If all elements are NaN
	  y = nan(class(x));
  end
else
  x = double(x); % Single failed in large arrays according to CAT12
  if nargin == 1              % Determine first nonsingleton dimension
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
			  half = floor(nnonnans(j)/2);
			  if 2*half == nnonnans(j)        % Average if even number of elements
				  y(:,j) = (x(half,j)+x(half+1,j))/2;
			  else
				  y(:,j) = x(half+1,j);
			  end
		  end
	  end
  else % There are no nans
      nCompare = s(dim);           % Number of elements used to generate a median
      half = floor(nCompare/2);    % Midway point, used for median calculation
      
      % If calculating along columns, use vectorized method with column
      % indexing.  Reshape at end to appropriate dimension.
      y = x(half+1,:);
      if 2*half == nCompare
          y = (x(half,:)+y)/2;
      end
      
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