function x_new = ps_scale(x, min_new, max_new)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

  x_min = min(x);
  x_new = x - x_min;
  x_max = max(x_new);
  x_new = x_new / x_max;
  x_new = x_new * (max_new - min_new) + min_new;

end

