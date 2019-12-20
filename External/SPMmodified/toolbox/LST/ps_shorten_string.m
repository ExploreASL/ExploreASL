function str_sh = ps_shorten_string(str, l)
%PS_SHORTEN_STR This function shortens the strint str to a desired length l
%
%   str_sh = ps_shorten_string(str, l) shortens string str to length l.
%
% Author: Paul Schmidt (www.statistical-modelling.de)

str_sh = str;
if numel(str) > l
    str_sh = ['...', str((end-l):end)];
end

end

