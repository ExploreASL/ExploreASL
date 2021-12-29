function T = xASL_stat_ticdf(P,nu)
% Student's t-distribution's inverse cumulative distribution function.
%
% FORMAT: T = xASL_stat_ticdf(P,nu)
%
% INPUT:
%   P        - Value of the CDF
%   nu       - Degrees of freedom of the CDF, either scalar or same size as T
%
% OUTPUT:
%   T        - Value of the inverse
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Calculates the inverse of cumulative distribution function of the Student's t-distribution for degrees of freedom nu at value P.
%
% EXAMPLE: T = xASL_stat_ticdf(0.1,5)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES: C. Walck. Hand-book on STATISTICAL DISTRIBUTIONS for experimentalists. 1996. University of Stockholm. 
%             http://inspirehep.net/record/1389910/files/suf9601.pdf
%             M. Abramowitz and I.A. Stegun. Handbook of mathematical Functions with Formulas, Graphs, and Mathematical Tables. 1964
%             http://people.math.sfu.ca/~cbm/aands/abramowitz_and_stegun.pdf
% __________________________________
% Copyright (C) 2015-2019 ExploreASL

% Admin
if isscalar(P) && (~isequal(size(P),size(nu)))
    P = P*ones(size(nu));
end

T = zeros(size(P));

% Set the trivial values
T(P<0) = NaN;
T(P>1) = NaN;
T(P==0) = -Inf;
T(P==1) = Inf;

% Use inversion of beta function for small nu
ii = find(nu < 1000);
if any(ii)
    Ploc = 2*P(ii)-1;
    nuloc = nu(ii);
	
	% For smaller values, calculate the inverse using the lower tail
    jj = (abs(Ploc)<0.5);
    X = zeros(size(Ploc));
    if any(jj)
        X(jj) = 1 - betaincinv(abs(Ploc(jj)),0.5,nuloc(jj)/2,'lower');
    end
    jj = ~jj;
    if any(jj)
        X(jj) = betaincinv(abs(Ploc(jj)),nuloc(jj)/2,0.5,'upper');
    end
    T(ii) = sign(Ploc).*sqrt(nuloc.*(1-X)./X);
end

% Abramowitz & Stegun formula 26.7.5 for large nu
ii = find(nu >= 1000);
if any(ii)
	Xp = -sqrt(2).*erfcinv(2*P(ii));
	g1 = (Xp.^3+Xp)/4;
	g2 = (5*Xp.^5+16*Xp.^3+3*Xp)/96;
	g3 = (3*Xp.^7+19*Xp.^5+17*Xp.^3-15*Xp)/384;
	g4 = (79*Xp.^9+776*Xp.^7+1482*Xp.^5-1920*Xp.^3-945*Xp)/92160;
	T(ii) = Xp + g1./nu(ii) + g2./(nu(ii).^2) + g3./(nu(ii).^3) + g4./(nu(ii).^4);
end
return
