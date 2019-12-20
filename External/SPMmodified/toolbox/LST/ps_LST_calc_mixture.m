function dens_csfgmwm = calc_mixture(Label, FLAIR, Lesion)

%     calc_mixture.m - Update parameters for gaussian mixture model, 
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

    CSF = FLAIR .* (Label < 1.5 & Lesion < 0.5);
    GM = FLAIR .* (Label >= 1.5 & Label < 2.5 & Lesion < 0.5);
    WM = FLAIR .* (Label >= 2.5 & Lesion < 0.5);
    
    sigma_csf = std(CSF(CSF > 0));
    sigma_gm = std(GM(GM > 0));
    sigma_wm = std(WM(WM > 0));
    
    mu_csf = mean(CSF(CSF > 0));
    mu_gm = mean(GM(GM > 0));
    mu_wm = mean(WM(WM > 0));
    
    n_csf = sum(CSF(:) > 0);
    n_gm = sum(GM(:) > 0);
    n_wm = sum(WM(:) > 0);

    PI_csf = n_csf / (n_csf + n_gm + n_wm);
    PI_gm = n_gm / (n_csf + n_gm + n_wm);
    PI_wm = n_wm / (n_csf + n_gm + n_wm);
    
    dens_csfgmwm = PI_csf .* ps_LST_dnorm(FLAIR, mu_csf, sigma_csf) + PI_gm .* ps_LST_dnorm(FLAIR, mu_gm, sigma_gm) + PI_wm .* ps_LST_dnorm(FLAIR, mu_wm, sigma_wm);
    %dens_csfgmwm(Lesion > 0.01) = 0;
    
return
