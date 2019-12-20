function [MU, SIGMA2, PI] = gmm_em(y, K)

%     gmm_em.m - Simple EM-algorithm for gaussian mixture model, 
%     part of the LST toolbox
%     Copyright (C) 2011  Paul Schmidt
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Initialization
PI = rand(1, K);
PI = PI ./ sum(PI);
MU = (max(y) - min(y)) .* rand(1, K);
SIGMA2 = repmat(var(y), 1, K);
stop_algorithm = 0;
ii = 0;
dens = zeros(numel(y), K);

while ~stop_algorithm
    ii = ii + 1;
    MU_old = MU;
    
    %% E-Step
    for kk = 1:K
        dens(:,kk) = PI(kk) .* ps_dnorm(y, MU(kk), sqrt(SIGMA2(kk)));
    end
    dens_sum = sum(dens, 2);
    P = dens ./ repmat(dens_sum, 1, K);

    %% M-step
    PI = mean(P);
    for kk = 1:K
        MU(kk) = sum(P(:,kk) .* y) ./ sum(P(:,kk));
        SIGMA2(kk) = var(y, P(:,kk)); 
    end

    %% check convergence
    if MU_old - MU < 0.0001
        stop_algorithm = 1;
    end

end

end