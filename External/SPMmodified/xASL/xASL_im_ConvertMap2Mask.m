function [ IMout ] = xASL_im_ConvertMap2Mask( IMin )
%xASL_im_ConvertMap2Mask Provides a robust way of conversion of 
% a continuous map to a binary mask, which can be used for lesions, ROIs,
% or tissue probability maps. Based on the assumption that a map should
% be thresholded at 50% to form a map, which is often the case for
% automatic segmentations

%% Explanation: taking the 0.5 could give a too strict mask, if there are a few high outliers
%  in case of interpolation overshoots or other causes for extreme values
%  (vascular peaks?). Masking at 0.5* the 0.95 percentile value instead of 1
%  (absolute max) helps, but doesn't work if a lot of low values have
%  occurred through smoothing from interpolation. A pragmatic solution is
%  too first make the strict mask, take the 0.95 percentile from this and
%  use half of its value.


% Example cases:
% Eg 1: in the case of a tissue probability map (SPM12/CAT12), this gives
% 0.5

% Eg 2: in the case of a downsampled (B-spline interpolated, native to
% standard space) this gives 0.5337 (simple half value would be 0.69)

% Eg 3: In case of already a mask, this gives 0.5 (i.e. leaves mask
% unchanged)

% Eg 4: In case of a probability segmentation (WMH PreDiva), this gives 0.5

% Eg 5: In case of a B-spline downsampled (native2MNI) probability
% segmentation (WMH PreDiva), this gives 0.5547 (simple half value would be
% 0.66)

% Eg 6: In case of the original functional ROI (fMRI activation), this
% gives 21.2 (simple half value would gives 29.8, real minimum of the ROI is 21.8))

% Eg 7: In case of 4rd degree B-spline interpolated functional ROI, this
% gives 21.2 (simple half value would give 29.6)

% Eg 8: In case of trilinear interpolated functional ROI, this
% gives 20.8 (simple half value would give 28.8)

% Eg 9: In case of nearest neighbor interpolated functional ROI, this
% gives 21.5 (simple half value would give 29.8)

% fprintf('%s\n','Map thresholded to mask, at "robust 50%" (i.e. 50% of the robust max (95%) ROI value');

if  xASL_stat_SumNan(IMin(:))>0

    %% Mask at 50%
    SimpleHalfValue                         = 0.5*max(IMin(:));
%     SimpleMask                              = IMin;
%     SimpleMask(SimpleMask<SimpleHalfValue)  = 0;

    % Exclude 5% highest values
%     MapValues       = sort(SimpleMask(SimpleMask>0)); % sorting all values within the map, assuming that >0 is within the mask
%     ValueM          = MapValues(round(0.95*length(MapValues))); % Clipping the highest 5% values
%     RobustHalfValue = 0.5*ValueM; % often, 0.5 is considered a safe threshold
%     IMout           = IMin>RobustHalfValue;

    % Bypassed the previous, to keep fair across all segmentations
    IMout           = IMin>SimpleHalfValue;
else
    IMout           = IMin;
end




end

