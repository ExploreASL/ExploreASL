function pth = ps_fullfile(varargin)
%ps_fullfile Like fullfile, but with the right file seperators.

pth = varargin{1};
for i = 2:nargin   
    pth = [pth, '/', varargin{i}];
end
pth = strrep(pth, '\', '/');

end

