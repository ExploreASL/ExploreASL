function [BLOB ROIvar] = Create3DROIblob( x, ASL, MaskChoice, PeakValue)
%Create3DROIblob Creates 3D Gaussian blob within ROI
% as artificial signal
% MaskChoice should be number corresponding to ROI (currently 1:14)

    %% Create BLOB
    MeanROI                 = xASL_stat_MeanNan(x.masks.Data.data(:,:,:,:,MaskChoice),4);
    BLOB                    = X_Y_Z_smoothing(single(MeanROI).*100000,8,8,8);
    BLOB                    = BLOB./max(BLOB(:)); % normalize
    BLOB                    = BLOB .* (BLOB>0.1); % confine to know "Ground Truth"

    %% Intensity
    BLOB                    = BLOB .* PeakValue; % set to PeakValue

    %% Get baseline variability throughout subjects
    ROImask                 = MeanROI > 0.8;

    for iSub= 1:size(ASL,4)
        temp                = squeeze(ASL(:,:,:,iSub));
        ROIvar(iSub,1)      = mean(mean(mean( temp(ROImask) )));
        clear temp
    end

    % Mean center Var
    ROIvar              = ROIvar./max(ROIvar(:));
    ROIvar              = (ROIvar+1-mean(ROIvar));

    % Pseudo-random shuffle
    ROIvar(:,2)         = rand(1,size(ROIvar,1));
    ROIvar              = sortrows(ROIvar,2);

end
