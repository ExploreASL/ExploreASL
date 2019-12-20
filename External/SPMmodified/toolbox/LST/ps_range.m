function r = ps_range(x)
%ps_LST_tlv   Check if the newest version of LST is installed.
%   Part of the LST toolbox, www.statistical-modeling.de/lst.html

r = diff([min(x); max(x)]);

end

