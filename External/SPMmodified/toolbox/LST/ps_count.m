function nx = count(x)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    x_sort = sort(x);
    x_sort_diff = diff(x_sort);
    x_unique = unique(x_sort);
    pos_diff = [0; find(x_sort_diff > 0); numel(x_sort)];
    nx = [x_unique, diff(pos_diff)]';

    %{ 
    %The old one
    labels = unique(x);
    nx = zeros(numel(labels), 1);
    for ii = 1:numel(labels)
        ff = find(x == labels(ii));
        nx(ii) = numel(ff);
        x(ff) = [];        
    end
    nx = [labels, nx]';
    %}
end

