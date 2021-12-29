function [b,CI,pval,stats] = xASL_stat_MultipleLinReg(X,Y,bIntercept)
% Calculates the multiple linear regression.
%
% FORMAT: [b,CI,pval,stats] = xASL_stat_MultipleLinReg(X,Y[,bIntercept])
%
% INPUT:
%   X          - A matrix NxK of K independent variables with N samples (REQUIRED)
%   Y          - A vector Nx1 of the dependent variable (REQUIRED)
%   bIntercept - Use intercept (OPTIONAL, DEFAULT = true);
%
% OUTPUT:
%   b        - Regression coefficients in a format [b1, b2, ...bk,a]
%   CI       - 95%-confidence interval for regression coefficients. Matrix of (K+1) x 2 values
%   pval     - P-values for each coefficient in the same order as for b. Vector of (K+1) x 1 values
%   stats    - A structure with additional results
%              stats.rSQ     - coefficient of determination R^2 describing the goodness of fit
%                            - R^2 of 0.95 indicates that 95% of the variation in Y is explained by the model
%              stats.rSQadj  - adjusted R^2
%              stats.F       - F-value
%              stats.pval    - pvalue of the global test of model adequacy
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: 
% Performs a multiple linear regression Y=b*X+a and provides the intercept and regression coefficients beta
% including their significance and confidence intervals. It calculates additionally the goodness of the fit.
%
% EXAMPLE: [b,CI,pval,stats] = xASL_stat_LinearRegress(X,Y)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES: C. Walck. Hand-book on STATISTICAL DISTRIBUTIONS for experimentalists. 1996. University of Stockholm. 
%             http://inspirehep.net/record/1389910/files/suf9601.pdf
%             http://home.iitk.ac.in/~shalab/regression/Chapter3-Regression-MultipleLinearRegressionModel.pdf
%             https://www.wessa.net/rwasp_multipleregression.wasp
% __________________________________
% Copyright (C) 2015-2019 ExploreASL
% 

% Admin
if nargin < 2 || isempty(X) || isempty(Y)
	error('Needs at least two input arguments.');
end

% Gets the number of samples
N = size(Y,1);

if size(Y,2) > 1
	error('Y needs to be Nx1');
end

if size(X,1) ~= N
	error('Y and X need to have the same number of samples.');
end

if nargin < 3 || isempty(bIntercept)
	bIntercept = true;
end

% Converts to double precision
X = double(X);
Y = double(Y);

% Reformulating the problem as Y = b*Xloc+a
if bIntercept
	Xloc = [X ones(N,1)];
else
	Xloc = X;
end

% And number of independent variables
K = size(Xloc,2);

if N-K < 2
	error('Too few samples for too much independent variables.');
end

% Solving in the ordinary least squares sense
%b = ((Xloc'*Xloc)^(-1))*(Xloc')*Y;
b = pinv(Xloc)*Y;

% The estimated values of y
Yest = Xloc*b;

% The residuals
epsilon = Y-Yest;

% Square of the residuals
epsilonSQ = epsilon'*epsilon;

% Set alpha for the confidence interval
alphaCI = 0.05;
sbSQ = epsilonSQ/(N-K)* ((Xloc'*Xloc)^(-1));

% Calculate the standard error
se = sqrt(diag(sbSQ));

% Confidence interval is
% bj - t(a/2,n-k)(sqrt(std^2Cjj))
t = xASL_stat_ticdf((1-alphaCI/2),(N-K));
CI = [b-t*se, b+t*se];


% Calculate the 2-tail p-values
% t = bi/se(i) ~ t(a/2,n-k-1)\
% b is k+1 because the intercept is added
v = N-K;
t = b./se;
pval = 2*xASL_stat_tcdf(t,v);

% Calculate the overall goodness of the fit
% R^2 = 1-e'e/(sum(yi-meany)^2)
stats.rSQ    = 1-epsilonSQ/sum((Y-mean(Y)).^2);
stats.rSQadj = 1-(1-stats.rSQ)*(N-1)/(N-K);
stats.F = ((sum((Y-mean(Y)).^2)-epsilonSQ)/(K-1)) / (epsilonSQ/(N-K));
stats.pval = xASL_stat_fcdf(1/stats.F,N-K,K-1); % Significance probability for regression
stats.se = se; % standard error

return
