function [H, P, W] = xASL_stat_ShapiroWilk(x, alpha)
% xASL_stat_ShapiroWilk performs the Shapiro-Wilk or Shapiro-Francia test of normality of a distribution.
%
% FORMAT: [H, P, W] = xASL_stat_ShapiroWilk(x[, alpha])
%
% INPUT:
%  X     Sample vector. (REQUIRED. Length at least 4)
%  ALPHA Significance level for the test (OPTIONAL. DEFAULT 0.05)
%
% OUTPUT:
%  H    For H=1 - reject the null hypothesis at the significance level ALPHA.%     H = 0 => Do not reject the null hypothesis at significance level ALPHA.
%       For H=0 - do not reject the null hypotehsis
%  P    P-value - the probability of observing X given that the null hypothesis is true = distribuition is normal. 
%  W    Normalized test statistics
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Performs the statistical test of normality - null hypothesis is that the sample is from normal
%              distribution with unspecified mean and variance. Based on the sample kurtosis it performs either 
%              Shapiro-Wilk (for platykurtic) or Shapiro-Francia test (for leptokurtic).
% EXAMPLE: xASL_stat_ShapiroWilk(x)
%          xASL_stat_ShapiroWilk(x, [], -1)
%          xASL_stat_ShapiroWilk(x, 0.001, -1)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% REFERENCES: References: Royston P. "Algorithm AS R94", Applied Statistics (1995) Vol. 44, No. 4.
%   AS R94
% __________________________________
% Copyright © 2015-2019 ExploreASL
%
% 2019-06-16 JP


%
%   The Shapiro-Wilk hypotheses are: 
%   Null Hypothesis:        X is normal with unspecified mean and variance.






% Admin
if nargin < 1 || isempty(x)
	error('xASL_stat_ShapiroWilk: Needs at least one input argument.');
end

% Remove NaNs
x = x(~isnan(x));
x = x(:);

if length(x) < 4
   error('xASL_stat_ShapiroWilk: At least four non-NaN samples are needed.');
end

% Checks alpha
if nargin < 2 || isempty(alpha)
	alpha = 0.05;
else
	if numel(alpha) > 1
		error('xASL_stat_ShapiroWilk: Alpha needs to be a scalar value.');
	end
	if alpha > 1 || alpha <= 0
		error('xASL_stat_ShapiroWilk: Alpha needs to be between 0 and 1.');
	end
end

x    = sort(x);
n    = length(x);
mi   = -sqrt(2) .* erfcinv(2*((1:n)'-0.375)/(n+0.25));
miSq = sum(mi.^2);

if xASL_stat_kurtosis(x) > 3
    % Kurtosis is higher than 3 -> Leptokurtic distribution
    % Perform Shapiro-Francia test

	weights = mi./sqrt(miSq);
	W       = (sum(weights.*x)^2) / sum((x-mean(x)).^2);
	
	nLog    = log(n);
	u1      = log(nLog)-nLog;
	u2      = log(nLog)+2/nLog;
	mu      = -1.2725 + 1.0521*u1;
    sigma   = 1.0308 - 0.26758*u2;
	
	% Normalized Shapiro-Francia statistics
	SF      = (log(1-W)-mu)/sigma;
	
	% P-value 
	P = 1 - (0.5*erfc(-SF./sqrt(2)));
	   
else
    % Kurtosis is higher than 3 -> Platykurtic distribution
    % Perform Shapiro-Wilke test
	u       = 1/sqrt(n);
	weights = zeros(n,1);
	
	weights(n) = -2.706056*(u^5) + 4.434685*(u^4) - 2.071190*(u^3) - 0.147981*(u^2) + 0.221157*u + mi(n)/sqrt(miSq);
    weights(1) = -weights(n);
	
	if n >= 6
		weights(n-1) = -3.582633*(u^5) + 5.682633*(u^4) - 1.752461*(u^3) - 0.293762*(u^2) + 0.042981*u + mi(n-1)/sqrt(miSq);
		weights(2)   = -weights(n-1);
		
		count = 3;
		eps = (miSq - 2*(mi(n)^2) - 2*(mi(n-1)^2)) / (1-2*(weights(n)^2)-2*(weights(n-1)^2));
	else
		count = 2;
		eps = (miSq - 2*(mi(n)^2)) / (1-2*(weights(n)^2));
	end
	
	weights(count:(n-count+1)) = mi(count:(n-count+1))/sqrt(eps);
	
	W = (sum(weights.*x)^2) / sum((x - mean(x)).^2);

	if (n <= 11)
		mu    = -0.0006714*(n^3) + 0.0250540*(n^2) - 0.39978*n + 0.54400;
		sigma = exp(-0.0020322*(n^3) + 0.0627670*(n^2) - 0.77857*n + 1.38220);
		gamma = 0.459*n - 2.273;
		
		SW    = -log(gamma-log(1-W));
	else
		nLog  = log(n);
		mu    = 0.0038915*(nLog^3) - 0.083751*(nLog^2) - 0.31082*nLog - 1.5861;
		sigma = exp( 0.0030302*(nLog^2) - 0.082676*nLog - 0.4803);
		
		SW    = log(1 - W);
	end

	% Normalize the Shapiro-Wilk statistic
    SW = (SW-mu)/sigma;
	
    % P-value
	P  = 1 - (0.5*erfc(-SW./sqrt(2)));
end

% Test the null hypothesis
H  = (alpha >= P);

end

function k = xASL_stat_kurtosis(x)
% Calculates kurtosis = E[(X-mu(X))^4] / (E[(X-mu(X))^2])^2

% Subtract mean from the 
xmu = x - xASL_stat_MeanNan(x);

k = xASL_stat_MeanNan(xmu.^4)./((xASL_stat_MeanNan(xmu.^2)).^2);

end