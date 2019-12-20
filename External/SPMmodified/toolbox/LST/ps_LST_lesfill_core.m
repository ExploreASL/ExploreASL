function filled = ps_LST_lesfill_core(img, prob, seg, indx_brain)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Correct seg, if image is not T1
t1 = 1; t2 = 0; f2 = 0;
%if ~((mean(img(seg == 1)) < mean(img(seg == 2))) && (mean(img(seg == 2)) < mean(img(seg == 3))))
if mean(img(prob > 0)) > mean(img(seg == 3))
    t1 = 0;    
    f2 = 1;
    y = img(seg == 1);
    csf_mean = mean(y);
    csf_sd = std(y);
    csf_ci = [csf_mean - 1.5 * csf_sd, csf_mean + 1.5 * csf_sd];
    seg_new = seg;
    seg_new(img > csf_ci(1) & img < csf_ci(2) & seg > 0) = 1;

    c_tmp = ps_count(seg_new(seg_new > 0));
    if any(c_tmp(2,:) ./ sum(c_tmp(2,:)) > .8)
        f2 = 0;
        t2 = 1;        
        seg_new = seg;
    end
    seg = seg_new;
end

% expand lesions if surrounding voxels are not CSF
% expand lesions if surrounding voxels are not CSF
binary = 1 .* (prob > 0);
csf = 1 .* (seg == 1);
indx_img = 0 .* img; indx_img(indx_brain) = indx_brain;
binary2 = binary;

indx_tmp = find(binary > 0);
nh_seg = getNeighborhood2(seg, indx_tmp, 1);
nh_indx = getNeighborhood2(indx_img, indx_tmp, 1);
tmp = nh_indx .* (nh_seg > 1);
binary(tmp(tmp > 0)) = 1;

% expand CSF
indx_tmp = find(seg == 1);
nh_seg = getNeighborhood2(seg, indx_tmp, 1);
nh_indx = getNeighborhood2(indx_img, indx_tmp, 1);
tmp = nh_indx .* (nh_seg > 1);
csf(tmp(tmp > 0)) = 1;

% Grow from the outside
indx_nowm = indx_brain(binary(indx_brain) > 0 | csf(indx_brain) > 0);    
filled = img; filled(indx_nowm) = 0;
indx_tmp = find(seg > 0 & filled > 0);
nh = getNeighborhood2(filled, indx_tmp, 1);
filled(indx_tmp(sum(nh > 0) < 3)) = 0;
binary(indx_tmp(sum(nh > 0) < 3)) = 1;

y = img(seg == 3);
gm_mean = mean(y);
gm_sd = std(y);
gm_ci = [gm_mean - 2 * gm_sd, gm_mean + 2 * gm_sd];

y = img(seg == 3 & filled > 0);
wm_mean = mean(y);
wm_sd = std(y);
wm_ci = [wm_mean - 2 * wm_sd, wm_mean + 2 * wm_sd];

y = img(seg == 1 & prob == 0);
csf_mean = mean(y);
csf_sd = std(y);
csf_ci = [csf_mean - 2 * csf_sd, csf_mean + 2 * csf_sd];


counter = 0; st = 0;
while ~st
    counter = counter + 1;
    indx_tmp = find(binary > 0);
    %numel(indx_tmp)
    if numel(indx_tmp) > 0
        nh = getNeighborhood2(filled, indx_tmp, 1);
        if t1
            nh(nh < csf_ci(2)) = 0;
        else
            if f2
                nh(nh < gm_ci(1) | nh > gm_ci(2)) = 0;
            else
                nh(nh < wm_ci(1) | nh > wm_ci(2)) = 0;
                %nh(nh > csf_ci(1)) = 0;
            end
        end
        nm = sum(nh) ./ sum(nh ~= 0);
        nm(isnan(nm)) = 0;
        if all(nh(:) == 0) | counter > 200
            st = 1;
        else
            filled(indx_tmp) = nm;% + normrnd(0, nanmean(nm) * .001, 1, numel(nm));
            %binary(indx_tmp(~isnan(nm))) = 0;
            binary(indx_tmp(sum(nh) ~= 0)) = 0;
        end
    else
        st = 1;
    end
end

if any(binary(:) > 0)
    csf(binary > 0) = 1;
end


indx_tmp = find(binary2 > 0);
nh = getNeighborhood2(filled, indx_tmp, 1);
tmp = bsxfun(@times, [filled(indx_tmp)'; nh], [1/3; repmat(1/9, 6, 1)]);
tmp(isnan(tmp)) = 0;
filled(indx_tmp) = sum(tmp);
filled(csf == 1 & binary2 == 0) = img(csf == 1 & binary2 == 0);

indx_tmp = find(seg > 0 & filled == 0);
if numel(indx_tmp) > 0
    st = 0;
    while ~st
        nh = getNeighborhood2(filled, indx_tmp, 1);
        filled(indx_tmp) = sum(nh) ./ sum(nh > 0);
        indx_tmp = find(seg > 0 & filled == 0);
        if isempty(indx_tmp)
            st = 1;
        end
    end    
end


end

