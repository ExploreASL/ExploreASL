function b = ps_bwlabeln(img)
%ps_bwlabeln   Label connected components in binary image
%   Part of the LST toolbox, www.statistical-modeling.de/lst.html
%
%   Sorry, there is no further documentation at this moment.
%

% Set borders of image to zero
img(:,:,1) = 0 .* img(:,:,1); img(:,:,end) = 0 .* img(:,:,end);
img(:,1,:) = 0 .* img(:,1,:); img(:,end,:) = 0 .* img(:,end,:);
img(1,:,:) = 0 .* img(1,:,:); img(end,:,:) = 0 .* img(end,:,:);

img = 1 .* (img ~= 0);
img2 = img;
nx = size(img, 1); ny = size(img, 2);
indx = find(img ~= 0);
indx1 = [indx - nx, indx - 1, indx + 1, indx + nx];
indx2 = indx - nx * ny;
indx3 = indx + nx * ny;
indx_nh = [indx1, indx2, indx3]';

indx_label = 0 .* indx;
indx_done = 0 .* indx;
label = 0;
lo = 1;
indx_label(lo) = label;
while any(indx_label == 0)    
    st = 0;
    label = label + 1;
    lo = find(indx_label == 0);
    lo = lo(1);
    indx_label(lo) = label;
    indx_done(lo) = 1;
    img2(indx(lo)) = 0;
    while ~st
        indx_nh_tmp = indx_nh(:,lo);
        [~, lo] = ismember(indx_nh_tmp(find(img2(indx_nh_tmp))), indx);
        lo = unique(lo);
        lo = lo(lo > 0);
        if isempty(lo)
            st = 1;
        else
            indx_label(lo) = label;
            indx_done(lo) = 1;
            img2(indx(lo)) = 0;
        end
    end 
end

b = img2;
b(indx) = indx_label;

end