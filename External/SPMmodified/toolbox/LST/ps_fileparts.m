function out = ps_fileparts(nam, type)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

type = sort(type);
[p, n, e] = fileparts(nam);
for i = 1:numel(type)
    if type == 1
        out = p;
    end
    if type == 2
        out = n;
    end
    if type == 3
        out = e;
    end
    if isequal(type, [1,2])
        out = fullfile(p, n);
    end
    if isequal(type, [2,3])
        out = [n, e];
    end    
end

end

