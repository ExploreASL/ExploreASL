Z_beta              =0.841621;              % Z(1-beta)  = Z(0.8)    = 0.84;
Z_alpha             =1.959964;              % Z(1-alpha) = Z(1-(0.05/2)) = Z(0.975) = 1.96;
% Z-score calculator -> www.fourmilab.ch/rpkp/experiments/analysis/zCalc.html

effect_size 					= 74.43-71.15 ;% (74.43-71.15)/74.43;

% effect_size                     = 9.9;
% variation_size_intrasubject     = 7.3;
VariationSize_betweensubject   = 14.8;

% (Z_beta+Z_alpha)^2 * (variation_size_intrasubject       /effect_size )^2
(Z_beta+Z_alpha)^2 * (VariationSize_betweensubject       /effect_size )^2