function q = ps_quantile(x, p)
%ps_LST_tlv   Check if the newest version of LST is installed.
%   Part of the LST toolbox, www.statistical-modeling.de/lst.html

sx = sort(x);
q = sx(ceil(numel(x) .* p));

end

