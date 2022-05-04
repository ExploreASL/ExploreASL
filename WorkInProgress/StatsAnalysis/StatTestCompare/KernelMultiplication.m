%% Single slice exercise for kernel multiplication to get difference


% 1 Create putamen BLOB

BLOB=BLOB(:,:,44); % select single slice
BLOB(BLOB==0)=NaN; % everything outside ROI = NaN
BLOB=BLOB./xASL_stat_MeanNan(BLOB(:)); % normalize to mean = 1

Trial = CBFImage{1}(:,:,44,1); % single data slice
Trial(isnan(BLOB))=NaN;     % select same putamen ROI
Trial=Trial./xASL_stat_MeanNan(Trial(:)); % normalize to 1

im1     = Trial;
im2     = Trial;

MultIm  = BLOB.*im2;
im3     = BLOB.*im2 ./ xASL_stat_MeanNan(MultIm(:)); % creates image with same mean, but new shape (BLOB-shape * data-shape)

xASL_stat_MeanNan(im1(:))-xASL_stat_MeanNan(im3(:))
im4=im1-im3;
xASL_stat_MeanNan(im4(:))
abs(xASL_stat_MeanNan(im1(:))-xASL_stat_MeanNan(im3(:)))/xASL_stat_MeanNan(im1(:))



