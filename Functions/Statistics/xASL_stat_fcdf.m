function p = xASL_stat_fcdf(F,M,N)
% F distribution function.
%
% FORMAT: F = xASL_stat_fcdf(F,M,N)
%
% INPUT:
%   F        - F-value for which the CDF is calculated
%   M,N      - Degrees of freedom of the CDF, either scalar or same size as T
%
% OUTPUT:
%   p        - Value of the CDF
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Calculates the cumulative distribution function of the F-distribution for degrees of freedom M,N at value F.
%
% EXAMPLE: p = xASL_stat_fcdf(F,M,N)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES: C. Walck. Hand-book on STATISTICAL DISTRIBUTIONS for experimentalists. 1996. University of Stockholm. 
%             http://inspirehep.net/record/1389910/files/suf9601.pdf
% __________________________________
% Copyright (C) 2015-2019 ExploreASL

% If degrees of freedom are given as a scalar value
if nargin < 3 || isempty(M) || isempty(N) || isempty(F)
	error('Need three input parameters.');
end
	
if isscalar(M)
	M = M*ones(size(F));
end

if isscalar(N)
	N = N*ones(size(F));
end

p = NaN(size(F));

ii = (F>0) & (M>0) & (N>0);
if sum(ii(:))
    jj = (ii & (F(ii) >= N(ii)./M(ii)));
	if sum(jj(:))
		p(jj) = betainc(N(jj)./(N(jj)+F(jj).*M(jj)), N(jj)/2, M(jj)/2, 'upper');
	end
	jj = (ii & (F(ii) < N(ii)./M(ii)));
    if sum(jj(:))
        p(jj) = betainc(M(jj).*F(jj)./(M(jj).*F(jj)+N(jj)), M(jj)/2, N(jj)/2, 'lower');
    end
end

p(F<=0) = 0;
return
