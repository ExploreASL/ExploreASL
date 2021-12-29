function F = xASL_stat_tcdf(T,nu)
% Student's t-distribution's cumulative distribution function.
%
% FORMAT: F = xASL_stat_tcdf(T,nu)
%
% INPUT:
%   T        - T-value for which the CDF is calculated
%   nu       - Degrees of freedom of the CDF, either scalar or same size as T
%
% OUTPUT:
%   F        - Value of the CDF
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Calculates the cumulative distribution function of the Student's t-distribution for degrees of freedom nu at value T.
%
% EXAMPLE: F = xASL_stat_tcdf(0.1,5)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES: C. Walck. Hand-book on STATISTICAL DISTRIBUTIONS for experimentalists. 1996. University of Stockholm. 
%             http://inspirehep.net/record/1389910/files/suf9601.pdf
% __________________________________
% Copyright (C) 2015-2019 ExploreASL

Tsq = T.^2;

% If degrees of freedom are given as a scalar value
if isscalar(nu)
	nu = nu*ones(size(T));
end

F = NaN(size(T));
% Calculate the CDF
ii = (nu < Tsq) & (nu>0) & (nu~=1) & (~isnan(T));
if any(ii(:))
	F(ii) = 0.5*betainc(nu(ii)./(nu(ii) + Tsq(ii)), nu(ii)/2, 0.5, 'lower');
end

ii = (nu>Tsq) & (nu>0) & (nu~=1) & (~isnan(T));
if any(ii(:))
	F(ii) = 0.5*betainc(Tsq(ii)./(nu(ii) + Tsq(ii)), 0.5, nu(ii)/2, 'upper');
end

F(T<0) = 1 - F(T<0);

ii = (nu==1) & (~isnan(T));
if any(ii(:))
	F(ii) = (T(ii)>0)+acot(-T(ii))/pi;
end

ii = (T==0) & (nu>0);
if any(ii(:))
	F(ii) = 0.5;
end

return