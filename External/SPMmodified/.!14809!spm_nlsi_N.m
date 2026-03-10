function [Ep,Eg,Cp,Cg,S,F,L] = spm_nlsi_N(M,U,Y)
% Bayesian inversion of a linear-nonlinear model of the form F(p)*G(g)'
% FORMAT [Ep,Eg,Cp,Cg,S,F,L]= spm_nlsi_N(M,U,Y)
%
% Generative model
%__________________________________________________________________________
% 
% M.IS - IS(p,M,U) A prediction generating function name; usually an 
%        integration scheme for state-space models of the form
%
%        M.f  - f(x,u,p,M) - state equation:  dxdt = f(x,u)
%
%        that returns hidden states - x; however, it can be any nonlinear
%        function of the inputs u. I.e., x = IS(p,M,U)
%
% M.G  - G(g,M) - linear observer: y = (x - M.x')*G(g,M)'
%
% M.FS - function name f(y,M) - feature selection
%        This [optional] function performs feature selection assuming the
%        generalized model y = FS(y,M) = FS(x*G',M) + X0*P0 + e
%
% M.x  - The expansion point for the states (i.e., the fixed point)
%
