function q = ps_qnorm(p, mu, sigma)
%ps_LST_qnorm   ---
%   Part of the LST toolbox, www.statistical-modeling.de/lst.html

q = mu + sigma * sqrt(2) .* erfinv(2 .* p - 1);

end

