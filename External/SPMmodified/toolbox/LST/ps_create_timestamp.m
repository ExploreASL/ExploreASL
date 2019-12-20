function ts = ps_create_timestamp
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
c = round(clock());
ts = num2str(c(1));
if c(2) < 10
    ts = [ts, '0', num2str(c(2))];
else
    ts = [ts, num2str(c(2))];
end
if c(3) < 10
    ts = [ts, '0', num2str(c(3))];
else
    ts = [ts, num2str(c(3))];
end
if c(4) < 10
    ts = [ts, '_0', num2str(c(4))];
else
    ts = [ts, '_', num2str(c(4))];
end
if c(5) < 10
    ts = [ts, '0', num2str(c(5))];
else
    ts = [ts, num2str(c(5))];
end
if c(6) < 10
    ts = [ts, '0', num2str(c(6))];
else
    ts = [ts, num2str(c(6))];
end

end

