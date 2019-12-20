function ps_LST_create_glass_brain(img, b, nam, or)
%ps_LST_create_glass_brain   Creates PNG for glass brains.
%   Part of the LST toolbox, www.statistical-modeling.de/lst.html
%
%   Sorry, there is no further documentation at this moment.
%
            
fl = or(end); 
or = or(1:3);
if any(or ~= [1 2 3])
    img = permute(img, or);    
    b = permute(b, or);      
end
if fl
    if exist('flip', 'builtin')
        img = flip(img, 2);
        b = flip(b, 2);
    else
        img = flipdim(img, 2);
        b = flipdim(b, 2);
    end
end
b = ps_set_border_zero(b);
indx_tmp = find(b > 0);
nh = getNeighborhood2(b, indx_tmp, 1);
border = 0 .* b;
border(indx_tmp(sum(nh == 0) > 0)) = 1;    

border_tmp = zeros(size(b, 2), size(b, 3));
border_tmp(:) = max(border, [], 1);
indx_tmp = find(border_tmp > 0);
nh = getNeighborhood2(border_tmp, indx_tmp, 0);
border1 = 0 .* border_tmp;
border1(indx_tmp(sum(nh == 0) > 0)) = 1;

border_tmp = zeros(size(b, 1), size(b, 3));
border_tmp(:) = max(border, [], 2);
indx_tmp = find(border_tmp > 0);
nh = getNeighborhood2(border_tmp, indx_tmp, 0);
border2 = 0 .* border_tmp;
border2(indx_tmp(sum(nh == 0) > 0)) = 1;

border_tmp = zeros(size(b, 1), size(b, 2));
border_tmp(:) = max(border, [], 3);
indx_tmp = find(border_tmp > 0);
nh = getNeighborhood2(border_tmp, indx_tmp, 0);
border3 = 0 .* border_tmp;
border3(indx_tmp(sum(nh == 0) > 0)) = 1;

ch11 = 0 .* border1; ch12 = 0 .* border2; ch13 = 0 .* border3; 
tmp = 1 .* (img > 0);
ch11(:) = mean(tmp, 1); ch12(:) = mean(tmp, 2); ch13(:) = mean(tmp, 3);
ch11 = (-1) .* ch11 + 1; ch12 = (-1) .* ch12 + 1; ch13 = (-1) .* ch13 + 1;
if sum(any(ch11 > 0 & ch11 < 1) > 0) > 0
    ch11(ch11 > 0 & ch11 < 1) = ps_scale(ch11(ch11 > 0 & ch11 < 1), 0, 1);
end
if sum(any(ch12 > 0 & ch12 < 1) > 0) > 0
    ch12(ch12 > 0 & ch12 < 1) = ps_scale(ch12(ch12 > 0 & ch12 < 1), 0, 1);
end
if sum(any(ch13 > 0 & ch13 < 1) > 0) > 0
    ch13(ch13 > 0 & ch13 < 1) = ps_scale(ch13(ch13 > 0 & ch13 < 1), 0, 1);
end
ch11(border1 > 0) = .5; ch12(border2 > 0) = .5; ch13(border3 > 0) = .5;

if exist('flip', 'builtin')
    imwrite(flip(ch11', 1), [nam, '_1.png'])
    imwrite(flip(ch12', 1), [nam, '_2.png'])
    imwrite(flip(ch13', 1), [nam, '_3.png'])
else
    imwrite(flipdim(ch11', 1), [nam, '_1.png'])
    imwrite(flipdim(ch12', 1), [nam, '_2.png'])
    imwrite(flipdim(ch13', 1), [nam, '_3.png'])
end
return

