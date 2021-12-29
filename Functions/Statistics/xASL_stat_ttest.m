function [H,P,CI,stats] = xASL_stat_ttest(X,M,alpha,tail,dim)
% One sample or paired-sample t-test of data from a normal distribution.
%
% FORMAT: [H,P,CI,stats] = xASL_stat_ttest(X[,M,alpha,tail,dim])
%
% INPUT:
%   X        - Input data - a vector or a N-dimensional matrix. NaNs are ignored
%   M        - [scalar] - mean value of the distribution according to the null hypothesis (OPTIONAL, DEFAULT 0)
%            - [matrix] of the same size as X - perfoms a t-test that the two matched samples (X and M) come
%              from two distributions with equal means
%   alpha    - Significance level of the test (OPTIONAL, DEFAULT 0.05)
%   tail     - a string (or a scalar) specifying if the test is single or two-tailed (OPTIONAL, DEFAULT 'both')
%            - 'both'  -  0 - two-tailed test
%            - 'right' -  1 - right-tailed test, testing if the mean > M
%            - 'left'  - -1 - left-tailed test, testing if the mean < M
%   dim      - The dimension along which to perform the test (OPTIONAL, DEFAULT - first non-singleton dimension)
%
% OUTPUT:
%   H        - Rejection of the null hypothesis
%            - H==1 means that the null hypothesis that the mean of the distribution of X is equal to zero can be rejected
%              at the given significance level Alpha
%            - H==0 the null hypothesis cannot be rejected at the given significance level
%   P        - P-value of the test
%   CI       - (1-alpha)-confidence interval for the mean of X
%   stats    - A structure with additional results
%            stats.tstat - t-statistics
%            stats.df    - degrees of freedom
%            stats.sd    - standard deviation of X
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Performs a t-test that the distribution of the input data X has a mean different from 0 (or from a
%              given mean M, or that the distributions X and Y have different means). A normal distribution of the data
%              with an unknown variance is assumed.
%
% EXAMPLE: H = xASL_stat_ttest(X)
%          [H,P,CI] = xASL_stat_ttest(X,0,0.01,'left',1)
%          [H,P,CI] = xASL_stat_ttest(X,Y,[]  ,'both',1)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES: C. Walck. Hand-book on STATISTICAL DISTRIBUTIONS for experimentalists. 1996. University of Stockholm. 
%             http://inspirehep.net/record/1389910/files/suf9601.pdf
% __________________________________
% Copyright (C) 2015-2019 ExploreASL

% Admin
if nargin < 2 || isempty(M)
    % Default mean is zero
    M = 0;
else
    % M is either a scalar or equally sized as X
    if isscalar(M) || isequal(size(M),size(X))
		X = X-M;
	else
        error('M has to be a scalar or equally sized as X');
    end
end

% Alpha is a scalar with values between 0 and 1
if nargin < 3 || isempty(alpha)
    alpha = 0.05;
elseif ~isscalar(alpha) || alpha <= 0 || alpha >= 1 || isnan(alpha)
    error('Significance level alpha needs to be a scalar between 0 and 1');
end

% By default, two-tailed test is performed
if nargin < 4 || isempty(tail)
    tail = 0;
else
	% Either tail is one of the given strings
	if ischar(tail)
		switch (tail)
			case 'left'
				tail = -1;
			case 'both'
				tail = 0;
			case 'right'
				tail = 1;
			otherwise
				tail = nan;
		end
	else
		% Or the correct string is given
		if ~isscalar(tail)
			tail = nan;
		end
	end
	
	% If outside of these values then report an error
	if ~sum(tail==[-1,0,1])
		error('Tail must be either -1, 0, or 1 or ''both'', ''right'', or ''left''');
	end
end

if nargin < 5 || isempty(dim)
    % Figure out which dimension mean will work along
	% By default the first non-singleton dimension
    dim = find(size(X) ~= 1, 1);
	if isempty(dim)
		dim = 1;
	end
end

% Identify the NaNs and calculate the sample size excluding NaNs
nans = isnan(X);
if any(nans(:))
    N = sum(~nans,dim);
else
    N = size(X,dim);
end

df = max(N-1,0);

% Calculate the population mean, standard deviation, and standard error of the mean
meanX = xASL_stat_MeanNan(X,dim);
stdX  = xASL_stat_StdNan(X,[],dim);
semX = stdX./sqrt(N);
tval = (meanX)./semX;

% Calculate the P-values
switch (tail)
	case  0
		% Two-tailed test
		P = 2 * xASL_stat_tcdf(abs(tval),df);
	case  1
		% Right-tailed test
		P = xASL_stat_tcdf(tval,df);
	case -1
		% Left-tailed test
		P = xASL_stat_tcdf(-tval,df);
end

% Calculate the confidence intervals
if nargout > 2
	switch (tail)
		case 0
			% Two-tailed test
			CIspread = xASL_stat_ticdf((1-alpha/2),df).*semX;
			CI = cat(dim,meanX-CIspread,meanX+CIspread);
		case 1
			% Right-tailed test
			CIspread = xASL_stat_ticdf(1-alpha,df).*semX;
			CI = cat(dim,meanX-CIspread,Inf(size(P)));
		case -1
			% Left-tailed test
			CIspread = xASL_stat_ticdf(1-alpha,df).*semX;
			CI = cat(dim,-Inf(size(P)),meanX+CIspread);
	end
end

% Create the outputs
% The rejection of the null hypothesis
H = double(P <= alpha);
H(isnan(P)) = nan;

% Add the mean to the confidence interval
if nargout > 2 && isscalar(M)
	CI = CI+M;
end

% Output the additional statisticas
if nargout > 3
	stats.tstat = tval;
	stats.sd    = stdX;
	stats.df    = df;

    if ~isequal(size(df),size(tval))
        stats.df = repmat(stats.df,size(tval));
    end
end

return
