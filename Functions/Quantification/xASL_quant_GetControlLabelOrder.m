function [ControlIm, LabelIm, OrderContLabl] = xASL_quant_GetControlLabelOrder(ASLTimeSeries)
%xASL_quant_GetControlLabelOrder Get control-label order of ASL time-series
%
% FORMAT: [ControlIm, LabelIm, OrderContLabl] = xASL_quant_GetControlLabelOrder(ASLTimeSeries)
%
% INPUT:
%   ASLTimeSeries - 4D matrix with multiple 3D ASL image volumes over time (REQUIRED)
%
% OUTPUT:
%   ControlIm     - 4D matrix with 3D ASL control image volumes
%   LabelIm       - 4D matrix with 3D ASL label image volumes
%   OrderContLabl - false for normal order (i.e. control-label)
%                   true for flipped order (i.e. label-control)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function automatically checks (and corrects if required)
%              the control and label order of ASL timeseries
%              based on the larger signal in control volumes.
%              It supposes that data is acquired in pairs. Works also for multiPLD but only for sequences
%              with alternative control/label or label/control order
%
% EXAMPLE:
% [~, ~, OrderContLabl] = xASL_quant_GetControlLabelOrder(ASLTimeSeries);
% [ControlIm, LabelIm] = xASL_quant_GetControlLabelOrder(ASLTimeSeries);
% 
% __________________________________
% Copyright (C) 2015-2021 ExploreASL

%% Get control-label order
ControlIm = ASLTimeSeries(:,:,:,1:2:end-1); % usual order
LabelIm = ASLTimeSeries(:,:,:,2:2:end-0);


% Check equality of n frames control & label
if size(ControlIm,4)~=size(LabelIm,4)
    error('Control and label image have unequal number of frames!');
end

%% Create time-average images
Ctrl = mean(ControlIm,4);
Lbl = mean(LabelIm,4);

%% Create brainmask
% First remove extreme values, we want to have average values
% within the brain
Ctrl = xASL_im_ClipExtremes(Ctrl, 0.99, 1, false); % 0.99 removes artifacts
% If this threshold doesn't work, it should be adjusted!

MaskIM = Ctrl>0.5.*max(Ctrl(:));

%% Masked images
Ctrl = Ctrl.*MaskIM;
Lbl = Lbl.*MaskIM;

%% Default order is [control label control label ...]
OrderContLabl = false;

TempPWI = Ctrl-Lbl; % [control label control label...]
% by default the ASL signal is positive in the control scan
% & negative in the label scan by the inversion
TempPWI(TempPWI==0) = NaN; % to enable taking median


% change control-label order if perfusion-weighted image is negative
% median is less sensitive to outliers
% it would be 0, if we hadn't changed zeros into NaNs
if xASL_stat_MedianNan(TempPWI(:))<0
    OrderContLabl = true;
    
    DummyIM = ControlIm;
    ControlIm = LabelIm;
    LabelIm = DummyIM;
end

end
