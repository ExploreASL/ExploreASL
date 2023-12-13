function [ControlIm, LabelIm, OrderContLabl, xQ] = xASL_quant_GetControlLabelOrder(ASLTimeSeries, xQ)
%xASL_quant_GetControlLabelOrder Get control-label order of ASL time-series
%
% FORMAT: [ControlIm, LabelIm, OrderContLabl, xQ] = xASL_quant_GetControlLabelOrder(ASLTimeSeries[, xQ])
%
% INPUT:
%   ASLTimeSeries - 4D matrix with multiple 3D ASL image volumes over time (REQUIRED)
%   xQ     - xQ field with input parameters containing the following subfields (OPTIONAL)
%                           - EchoTime: TE vector for unsubtracted images
%                           - Initial_PLD: PLD vector for unsubtracted images
%                           - LabelingDuration: : LD vector for unsubtracted images
%
% OUTPUT:
%   ControlIm     - 4D matrix with 3D ASL control image volumes
%   LabelIm       - 4D matrix with 3D ASL label image volumes
%   OrderContLabl - false for normal order (i.e. control-label)
%                   true for flipped order (i.e. label-control)
%   xQ             - xQ field with several output parameters - only calculated if xQ is on input
%                           - EchoTime_PWI4D: TE vector after subtraction
%                           - InitialPLD_PWI4D: PLD vector after subtraction
%                           - LabelingDuration_PWI4D: LD vector after subtraction
%                           - EchoTime_Control4D: TE vector for all control images
%                           - InitialPLD_Control4D: PLD vector for all control images

% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function automatically checks (and corrects if required) the control and label order of ASL timeseries
%              based on the larger signal in control volumes. It supposes that data is acquired in pairs. Works also for multiPLD but only for sequences
%              with alternative control/label or label/control order. It also calculated PLD, Labeling Duration and EchoTime vector for the output
%
% EXAMPLE:
% [~, ~, OrderContLabl] = xASL_quant_GetControlLabelOrder(ASLTimeSeries);
% [ControlIm, LabelIm, xQ] = xASL_quant_GetControlLabelOrder(ASLTimeSeries, xQ);
% 
% __________________________________
% Copyright (C) 2015-2023 ExploreASL

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

if nargin >= 2 && ~isempty(xQ)
	% xQ is provided, calculate the output parameters
	% Skip every other value in the vectors as they were stored for both control and label images
	xQ.EchoTime_PWI4D = xQ.EchoTime(1:2:end);
	xQ.InitialPLD_PWI4D = xQ.Initial_PLD(1:2:end);
	xQ.LabelingDuration_PWI4D = xQ.LabelingDuration(1:2:end);
	% In this case, the vectors are equal for PWI4D and Control4D
	xQ.EchoTime_Control4D = xQ.EchoTime_PWI4D;
	xQ.InitialPLD_Control4D = xQ.InitialPLD_PWI4D;
	xQ.LabelingDuration_Control4D = xQ.LabelingDuration_PWI4D;
end

end
