function y = prctile_fct(varargin)
%PRCTILE Percentiles of a sample.
%   Y = PRCTILE(X,P) returns percentiles of the values in X.  P is a scalar
%   or a vector of percent values.  When X is a vector, Y is the same size
%   as P, and Y(i) contains the P(i)-th percentile.  When X is a matrix,
%   the i-th row of Y contains the P(i)-th percentiles of each column of X.
%   For N-D arrays, PRCTILE operates along the first non-singleton
%   dimension.
%
%   Y = PRCTILE(X,P,'all') calculates percentiles of all the elements in X.
%   The smallest dimension index of Y has length LENGTH(P)
%
%   Y = PRCTILE(X,P,DIM) calculates percentiles along dimension DIM. The
%   DIM'th dimension of Y has length LENGTH(P).
%
%   Y = PRCTILE(X,P,VECDIM) calculates percentiles of elements of X based
%   on the dimensions specified in the vector VECDIM. The smallest dimension
%   index specified in VECDIM has length LENGTH(P).
%
%   Y = PRCTILE(...,'PARAM1',val1,'PARAM2',val2,...) specifies optional
%   parameter name/value pairs:
%
%   'Method'  - 'exact' (default) to compute by sorting as explained below. 
%               'approximate' to use an approximation algorithm based on 
%               t-digests.
%
%   Percentiles are specified using percentages, from 0 to 100.  For an N
%   element vector X, PRCTILE computes percentiles as follows:
%      1) The sorted values in X are taken as the 100*(0.5/N), 100*(1.5/N),
%         ..., 100*((N-0.5)/N) percentiles.
%      2) Linear interpolation is used to compute percentiles for percent
%         values between 100*(0.5/N) and 100*((N-0.5)/N)
%      3) The minimum or maximum values in X are assigned to percentiles
%         for percent values outside that range.
%
%   PRCTILE treats NaNs as missing values, and removes them.
%
%   Examples:
%      y = prctile(x,50); % the median of x
%      y = prctile(x,[2.5 25 50 75 97.5]); % a useful summary of x
%
%   See also IQR, MEDIAN, NANMEDIAN, QUANTILE.

%   Copyright 1993-2018 The MathWorks, Inc.

par = inputParser();
par.addRequired('x');
par.addRequired('p');
par.addOptional('dim',1,@(x) isnumeric(x) || validateDimAll(x));
par.addParameter('Method', 'exact');
par.addParameter('Delta',1e3);
par.addParameter('RandStream',[]);

par.parse(varargin{:});

x = par.Results.x;
p = par.Results.p;
dim = par.Results.dim;
method = par.Results.Method;
delta = par.Results.Delta;
rs = par.Results.RandStream;

if ~isvector(p) || numel(p) == 0 || any(p < 0 | p > 100) || ~isreal(p)
    error(message('stats:prctile:BadPercents'));
end
%validatestring(method,["exact", "approximate"]);
validatestring(method,{'exact', 'approximate'});

% Make sure we are working in floating point to avoid rounding errors.
if isfloat(x)
    castOutput = false;
elseif isinteger(x)
    % integer types are up-cast to either double or single and the result
    % is down-cast back to the input type
    castOutput = true;
    outType = internal.stats.typeof(x);
%    if ismember(outType, ["int8" "uint8" "int16" "uint16"])
    if ismember('int8', {'int8' 'uint8' 'int16' 'uint16'})
        % single precision is enough
        x = single(x);
    else
        % Needs double precision
        x = double(x);
    end
else
    % All other types (e.g. char, logical) are cast to double and the result is
    % double
    castOutput = false;
    x = double(x);
end

% Figure out which dimension prctile will work along.
sz = size(x);
if ismember('dim', par.UsingDefaults) 
    dim = find(sz ~= 1,1);
    if isempty(dim)
        dim = 1; 
    end
    dimArgGiven = false;
else
    if validateDimAll(dim)
        dim = 1:ndims(x);
    end
    % Checking validity of dim
    if isnumeric(dim)
        if ~isreal(dim) || any(floor(dim) ~= ceil(dim)) || any(dim < 1) || any(~isfinite(dim))
            error(message('MATLAB:getdimarg:invalidDim'));
        end
        if ~isscalar(dim) && ~all(diff(sort(dim)))
            error(message('MATLAB:getdimarg:vecDimsMustBeUniquePositiveIntegers'));
        end
    else
        error(message('MATLAB:getdimarg:invalidDim'));
    end  
    % Permute the array so that the requested dimension is the first dim.
    nDimsX = ndims(x);
    dim = sort(dim);
    perm = [dim setdiff(1:max(nDimsX,max(dim)),dim)];
    x = permute(x, perm);
    sz = size(x);
    dimArgGiven = true;
end

% If X is empty, return all NaNs.
if isempty(x)
    if isequal(x,[]) && ~dimArgGiven
        y = nan(size(p),'like',x);
    else
        if dimArgGiven
            work_dim = 1:numel(dim);
        else
            work_dim = dim;
        end
        szout = sz; 
        szout(work_dim) = 1;
        szout(work_dim(1)) = numel(p);
        y = nan(szout,'like',x);
    end

else
    % Drop X's leading singleton dims, and combine its trailing dims.  This
    % leaves a matrix, and we can work along columns.
    if dimArgGiven
        work_dim = 1:numel(dim);
    else
        work_dim = dim;
    end
    work_dim = work_dim(work_dim <= numel(sz));
    nrows = prod(sz(work_dim));
    ncols = numel(x) ./ nrows;
    x = reshape(x, nrows, ncols);
    
    if strcmpi(method,'exact')
        x = sort(x,1);
        n = sum(~isnan(x), 1); % Number of non-NaN values in each column
        
        % For columns with no valid data, set n=1 to get nan in the result
        n(n==0) = 1;
        
        % If the number of non-nans in each column is the same, do all cols at once.
        if all(n == n(1))
            n = n(1);
            if isequal(p,50) % make the median fast
                if rem(n,2) % n is odd
                    y = x((n+1)/2,:);
                else        % n is even
                    y = (x(n/2,:) + x(n/2+1,:))/2;
                end
            else
                y = interpColsSame(x,p,n);
            end
            
        else
            % Get percentiles of the non-NaN values in each column.
            y = interpColsDiffer(x,p,n);
        end
    else
        td = TDigestArray(x, delta, rs);
        y = internal.stats.bigdata.tdigestICDF(td,p,min(x,[],1),max(x,[],1));
    end

    % Reshape Y to conform to X's original shape and size.
    szout = sz; 
    szout(work_dim) = 1;
    szout(work_dim(1)) = numel(p);
    y = reshape(y,szout);
end
% undo the DIM permutation
if dimArgGiven
     y = ipermute(y,perm);  
end

% If X is a vector, the shape of Y should follow that of P, unless an
% explicit DIM arg was given.
if ~dimArgGiven && isvector(x)
    y = reshape(y,size(p)); 
end

if castOutput
    y = cast(y, outType);
end


function y = interpColsSame(x, p, n)
%INTERPCOLSSAME An aternative approach of 1-D linear interpolation which is
%   faster than using INTERP1Q and can deal with invalid data so long as
%   all columns have the same number of valid entries (scalar n).

% Make p a column vector. Note that n is assumed to be scalar.
if isrow(p)
    p = p';
end

% Form the vector of index values (numel(p) x 1)
r = (p/100)*n;
k = floor(r+0.5); % K gives the index for the row just before r
kp1 = k + 1;      % K+1 gives the index for the row just after r
r = r - k;        % R is the ratio between the K and K+1 rows

% Find indices that are out of the range 1 to n and cap them
k(k<1 | isnan(k)) = 1;
kp1 = bsxfun( @min, kp1, n );

% Use simple linear interpolation for the valid percentages
y = (0.5+r).*x(kp1,:)+(0.5-r).*x(k,:);

% Make sure that values we hit exactly are copied rather than interpolated
exact = (r==-0.5);
if any(exact)
    y(exact,:) = x(k(exact),:);
end

% Make sure that identical values are copied rather than interpolated
same = (x(k,:)==x(kp1,:));
if any(same(:))
    x = x(k,:); % expand x
    y(same) = x(same);
end

function y = interpColsDiffer(x, p, n)
%INTERPCOLSDIFFER A simple 1-D linear interpolation of columns that can
%deal with columns with differing numbers of valid entries (vector n).

[nrows, ncols] = size(x);

% Make p a column vector. n is already a row vector with ncols columns.
if isrow(p)
    p = p';
end

% Form the grid of index values (numel(p) x numel(n))
r = (p/100)*n;
k = floor(r+0.5); % K gives the index for the row just before r
kp1 = k + 1;      % K+1 gives the index for the row just after r
r = r - k;        % R is the ratio between the K and K+1 rows

% Find indices that are out of the range 1 to n and cap them
k(k<1 | isnan(k)) = 1;
kp1 = bsxfun( @min, kp1, n );

% Convert K and Kp1 into linear indices
offset = nrows*(0:ncols-1);
k = bsxfun( @plus, k, offset );
kp1 = bsxfun( @plus, kp1, offset );

% Use simple linear interpolation for the valid percentages.
% Note that NaNs in r produce NaN rows.
y = (0.5-r).*x(k) + (0.5+r).*x(kp1);

% Make sure that values we hit exactly are copied rather than interpolated
exact = (r==-0.5);
if any(exact(:))
    y(exact) = x(k(exact));
end

% Make sure that identical values are copied rather than interpolated
same = (x(k)==x(kp1));
if any(same(:))
    x = x(k); % expand x
    y(same) = x(same);
end

function bool = validateDimAll(dim)
bool = ((ischar(dim) && isrow(dim)) || ...
     (isstring(dim) && isscalar(dim) && (strlength(dim) > 0))) && ...
     strncmpi(dim,'all',max(strlength(dim), 1));
