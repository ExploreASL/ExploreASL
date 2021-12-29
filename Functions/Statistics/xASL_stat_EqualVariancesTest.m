function [resTest, P] = xASL_stat_EqualVariancesTest(X, alpha, type)
% Test for the equality of variances.
%
% FORMAT: [resTest, P] = xASL_stat_EqualVariancesTest(X[, alpha, type])
%
% INPUT:
%   X        - Input data (Nx2 matrix, first column are the samples, second column group identificators 1,2,3...)
%   alpha    - Significance level (DEFAULT 0.05)
%   type     - 'BrownForsythe' (DEFAULT) for the Brown-Forsythe test
%              'Levene' for the Levene's test
% OUTPUT:
%   resTest  - Result of the test (1 if any of the variances are not equal)
%   P        - P-value of the F-test
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Brown-Forsythe or Levene's test for equality of variances. The response variable is
%              transformed (yij = abs(xij - median(xj)) for Brown-Forsythe and yij = abs(xij - mean(xj))
%              for Levene's test). And then runs a one-way ANOVA F-test to check if the variances are equal.
%
% EXAMPLE: X = rand(10000,1);
%          X(1:5000,2) = 1;
%          X(5001:end,2) = 2;
%          [res,p] = xASL_stat_EqualVariancesTest(X,0.05,'Levene')
%          [res,p] = xASL_stat_EqualVariancesTest(X,[],'BrownForsythe')
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES: Brown, M. B. and Forsythe, A. B. (1974), Robust Tests for the Equality of Variances. 
%                 Journal of the American Statistical Association, 69:364-367.
%             Zar, J. H. (1999), Biostatistical Analysis (2nd ed.). NJ: Prentice-Hall, Englewood Cliffs. p. 180. 
% __________________________________
% Copyright (C) 2015-2019 ExploreASL
%
% 2019-04-24 JP, HM

% Admin
if nargin < 2 || isempty(alpha)
   alpha = 0.05;
end 

if nargin < 3 || isempty(type)
	type = 'BrownForsythe';
end

if ~strcmp(type,'BrownForsythe') && ~strcmp(type,'Levene')
	error('xASL_stat_EqualVariancesTest: Type needs to be BrownForsythe or Levene.');
end

% Number of groups
k=max(X(:,2));

% Transformed response variable Yij= |Xij-E(Xj)|
Y=X;
for i=1:k
   ind = (Y(:,2) == i);
   if strcmp(type,'BrownForsythe')
	   Y(ind,1) = abs(Y(ind,1) - median(Y(ind,1)));
   else
	   Y(ind,1) = abs(Y(ind,1) - mean(Y(ind,1)));
   end
end

%Analysis of variance procedure.
NX = zeros(1,k);
YGrpMean = zeros(1,k);
YTotMean = mean(Y(:,1));
YCor = zeros(size(Y,1),1);
for i=1:k
   ind = Y(:,2)==i;
   NX(1,i) = length(Y(ind,1));
   YGrpMean(1,i) = mean(Y(ind,1));
   YCor(ind,1) = Y(ind,1) - YGrpMean(1,i);
end

d1=k-1; %sample degrees of freedom.
d2=length(Y(:,1))-k; %error degrees of freedom.
SSA = sum(NX.*((YGrpMean-YTotMean).^2));
SSE = sum(YCor.^2);
F = d2/d1*SSA/SSE;

% F follows the F-distribution
P = 1 - fCumDistrFun(F,d1,d2);  

resTest= (P<alpha);

end

function p = fCumDistrFun(x,d1,d2)
% F-cumulative distribution function for finite d1,d2 > 0. For infinite x, returns 1.

% Initialize P to zero.
p = zeros(size(x));

% Return nans for badly defined parameters
if isnan(d1) || isnan(d2) || d1<=0 || d2<= 0 || (~isfinite(d1)) || (~isfinite(d2))
	p(:) = NaN;
	return;
end

% Catch Inf and NaNs and below zeros
p(isnan(x)) = NaN;
p(x==Inf) = 1;
p(x<0) = NaN;

% Calculate only for those
ind1 = (x > 0 & isfinite(x) & ~isnan(x));
if any(ind1)
	ind2 = d2 > x(ind1)*d1;
	
	% I_(d1x/(d1x+d2)) (d1/2,d2/2)
	if any(ind2)
		indLoc = ind1(ind2);
        p(indLoc) = betainc((d1*x(indLoc))./(d1*x(indLoc)+d2), d1/2, d2/2,'lower');
	end

	% I_(d2/(d2+d1x)) (d2/2,d1/2)
	if any(~ind2)
		indLoc = ind1(~ind2);
        p(indLoc) = betainc(d2./(d2+x(indLoc).*d1), d2/2, d1/2,'upper');
	end
end
 
end
