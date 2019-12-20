function img = ps_set_border_zero(img)

img(:,:,1) = 0 .* img(:,:,1); img(:,:,end) = 0 .* img(:,:,end);
img(:,1,:) = 0 .* img(:,1,:); img(:,end,:) = 0 .* img(:,end,:);
img(1,:,:) = 0 .* img(1,:,:); img(end,:,:) = 0 .* img(end,:,:);

end

