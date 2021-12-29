function [H,P,CI,stats] = xASL_stat_ttest2(X,Y,alpha,tail,vartype,dim)
% Two sample unpaired t-test of data from a normal distribution.
%
% FORMAT: [H,P,CI,stats] = xASL_stat_ttest2(X,Y[,alpha,tail,vartype,dim])
%
% INPUT:
%   X        - Input data - a vector or a N-dimensional matrix. NaNs are ignored
%   Y        - Second input data - a vector or a N-dimensional matrix. NaNs are ignored. 
%              Size must be similar to X. Except for the working dimension - this can differ
%   alpha    - Significance level of the test (OPTIONAL, DEFAULT 0.05)
%   tail     - a string (or a scalar) specifying if the test is single or two-tailed (OPTIONAL, DEFAULT 'both')
%            - 'both'  -  0 - two-tailed test
%            - 'right' -  1 - right-tailed test, testing if the mean > M
%            - 'left'  - -1 - left-tailed test, testing if the mean < M
%   vartype  - A string describing the variances of X a Y (OPTIONAL, DEFAULT 'equal')
%            - 'equal' - assumes equal variances of X and Y
%            - 'unequal' - assumes unequal variances of X and Y
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
%            stats.sd    - standard deviation of X. A vector is given for the unequal variances
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Performs a unpaired t-test that the distribution of the input data X has a mean different from that of Y. 
%              A normal distribution of the data with an unknown variance is assumed.
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
if nargin < 2
    error('Two samples need to be given.');
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

% Equal variances are assumed by default
if nargin < 5 || isempty(vartype)
    vartype = 1;
else
	if ischar(vartype) && (size(vartype,1)==1)
		switch(vartype)
			case 'equal'
				vartype = 1;
			case 'unequal'
				vartype = 2;
			otherwise
				error('vartype must be one of the strings ''equal'' or ''unequal''.');
		end
	else
		error('vartype must be one of the strings ''equal'' or ''unequal''.');
	end
end

if nargin < 6 || isempty(dim)
    % By default, use the first non-singleton dimension
	% First try X
	dim = find(size(X) ~= 1, 1);
	% Then try Y
	if isempty(dim)
		dim = find(size(Y) ~= 1, 1);
	end
	% Then use the first dimension
	if isempty(dim)
		dim = 1;
	end
end

% The data must have the same dimension, except for the working dimension - there, a different number of samples can be given
NX = size(X); 
NY = size(Y);
NY(dim) = NX(dim);
if ~isequal(NX,NY)
    error('Dimensions of X and Y must match in non-working dimensions.');
end

% Identify the NaNs and calculate the sample size excluding NaNs
nansX = isnan(X);
if any(nansX(:))
    NX = sum(~nansX,dim);
else
    NX = size(X,dim);
end
nansY = isnan(Y);
if any(nansY(:))
    NY = sum(~nansY,dim);
else
	NY = size(Y,dim); 
end

varX    = xASL_stat_VarNan(X,[],dim);
varY    = xASL_stat_VarNan(Y,[],dim);
dmeanXY = xASL_stat_MeanNan(X,dim) - xASL_stat_MeanNan(Y,dim);

switch (vartype)
	% Equal variances
	case 1 
		% Degrees of freedom
		df = NX + NY - 2;
		% Pooled variance
		stdXY = sqrt(((NX-1).*varX + (NY-1).*varY)./df);
		semXY = stdXY.*sqrt(1./NX+1./NY);
		tval  = dmeanXY./semXY;
	
	% Unequal variances
	case 2 
		semXY = sqrt(varX./NX + varY./NY);
		tval = dmeanXY./semXY;
		df = (varX./NX + varY./NY) .^2 ./ ((varX./NX).^2 ./(NX-1) + (varY./NY).^2 ./(NY-1));
		if semXY == 0
			df = 1; 
		end
end

% Calculate the P-values
switch (tail)
	case 0
		P = 2*xASL_stat_tcdf(abs(tval),df);
	case 1
		P = xASL_stat_tcdf(tval,df);
	case -1
		P = xASL_stat_tcdf(-tval,df);
end

if nargout > 2
	switch (tail)
		case 0
			CIspread = xASL_stat_ticdf(1-alpha./2,df).*semXY;
			CI = cat(dim, dmeanXY-CIspread, dmeanXY+CIspread);
		case 1
			CIspread = xASL_stat_ticdf(1-alpha, df).*semXY;
			CI = cat(dim, dmeanXY-CIspread,Inf(size(P)));
		case -1
			CIspread = xASL_stat_ticdf(1 - alpha, df).*semXY;
			CI = cat(dim, -Inf(size(P)),dmeanXY+CIspread);
	end
end

% Create the outputs
% The rejection of the null hypothesis
H = double(P <= alpha);
H(isnan(P)) = nan;

% Output the additional statisticas
if nargout > 3
	stats.tstat = tval;
	if vartype == 1
		stats.sd = stdXY;
	else
		stats.sd = sqrt(cat(dim, varX, varY));
	end
	stats.df    = df;

    if ~isequal(size(df),size(tval))
        stats.df = repmat(stats.df,size(tval));
    end
end

return
