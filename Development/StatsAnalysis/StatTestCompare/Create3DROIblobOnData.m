function [BLOB CBFout CBF1 GoldPos] = Create3DROIblobOnData( x, CBFin, CBF1, MaskChoice, ChangeFactor, ShapeFactor, ActivationVolume)
%Create3DROIblobOnData Imposes 4D Gaussian blob on ROI, for consistent change (artificial signal)
% This script takes the original signal distribution, multiplies it with an ROI convolved with 
% an 8mm isotropic Gaussian kernel to have a mean signal change (ChangeFactor) that follows the
% original data distribution but also the ROI shape (largest in center/skeleton and decreasing towards the ROI edges)
% MaskChoice should be number corresponding to ROI (currently 1:14)
% CBF1 = original group

% ROI = GM	WM	L-ICA	R-ICA	POS	6 Caudate	Cerebellum	8 Frontal	Insula	10 Occipital	Parietal Putamen	Temporal 14	Thalamus

    % Default: ActivationVolume = 1 (activation spread over 100% of ROI); ShapeFactor = 1/3 (<0 = platykurtic, >1 = leptokurtic); ChangeFactor = 1.05 (5 % effect size)
    


    %% Create BLOB
    MeanROI                 = xASL_stat_MeanNan(x.masks.Data.data(:,:,:,:,MaskChoice),4);
    BLOB                    = X_Y_Z_smoothing(single(MeanROI).*100000,8,8,8);

    %% Here, the Gaussian ROI-shaped cone is iteratively eroded to end up with a percentage of initial volume
    ThrFactor               = 0.01; % initial thresholding (erosion) factor
    GoldPos                 = BLOB > (ThrFactor*max(BLOB(:)));          % Ground Truth
    BLOB(~GoldPos)          = NaN;                                      % restrict BLOB/cone to Ground Truth

    % current volume
    TotalVol                = xASL_stat_SumNan(GoldPos(:));
    ActivVol                = round(ActivationVolume .* TotalVol);      % volume we wish to obtain
    CurrentVol              = sum(sum(sum(~isnan(BLOB))));

    % the loop
    while CurrentVol>ActivVol
        ThrFactor           = ThrFactor+0.01;
        GoldPos             = BLOB > (ThrFactor*max(BLOB(:)));          % reset Ground Truth
        BLOB(~GoldPos)      = NaN;                                      % restrict BLOB/cone to Ground Truth
        CurrentVol          = sum(sum(sum(~isnan(BLOB))));
    end
    
    %% Shape (lepto-kurtic vs. platy-kurtic)
    % Here without quantification of the shape, just taking powering it with factor that makes more flat/platykurtic (0.001) or more lepto-kurtic (10)
    BLOB                    = BLOB.^ShapeFactor;

    %% Normalize intensity within the new shape to 1, and multiply by requested "change"/effect size
    BLOB                    = BLOB./xASL_stat_MeanNan(BLOB(:));               % normalize to mean of 1
    BLOB                    = repmat(BLOB,[1 1 1 size(CBFin,4)]);   % repeat for all subjects
    Masks                   = repmat(GoldPos,[1 1 1 size(CBFin,4)]);

    CBFin(~Masks)           = NaN;
    CBF1(~Masks)            = NaN;

%     if      sum(sum(sum(sum( ~isnan(BLOB) )))) ~=  sum(sum(sum(sum( ~isnan(CBFin) )))) -> error could be related to existing NaNs
%             error('BLOB & masked data don"t have same voxel count!!!)');
%     end
    
    MultIm                  = BLOB.*CBFin;
    CBFout                  = (BLOB.*CBFin) ./ xASL_stat_MeanNan(MultIm(:)) .* xASL_stat_MeanNan(CBFin(:)) .* ChangeFactor;

    BLOB                    = BLOB(:,:,:,1); % restricts memory usage
    
   

%% Old visualization tools
%  
% BLOB(isnan(BLOB))           = 0;
% dip_image(50.*background_mask(:,:,44)+BLOB(:,:,44,1))
% BLOB(BLOB==0)               = NaN;
% Plot2Dsurface(CropNaN( BLOB(:,:,44,1), 0.1 ))
% TempData(isnan(TempData))   = 0;
% dip_image(50.*background_mask(:,:,44)+TempData(:,:,44,1))
% Plot2Dsurface(CropNaN( TempData(:,:,44,1), 0.1 ))
% CBFout(isnan(CBFout))       = 0;
% dip_image(50.*background_mask(:,:,44)+CBFout(:,:,44,1))
% Plot2Dsurface(CropNaN( CBFout(:,:,44,1), 0.1 ))
% 
% Plot2Dsurface(BLOB(:,:,44,1))
% Plot2Dsurface(GoldPos(:,:,44,1))
% dip_image(GoldPos(:,:,44,1))

%     Plot2Dsurface(TempData(:,:,44,1))
%     Plot2Dsurface(CBFout(:,:,44,1))
%     dip_image(GoldPos(:,:,44,1))    




end
