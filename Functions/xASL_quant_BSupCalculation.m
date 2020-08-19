function SignalPercentage = xASL_quant_BSupCalculation(BSupTime, PresaturationTime, T1Time, SliceTime, PathGraph)
% xASL_quant_BSupCalculation calculates the percentage of the remaining signal after a Bsup correction
%
% FORMAT: signalPercentage = xASL_quant_BSupCalculation(BSupTime[, PresaturationTime, T1Time, SliceTime, PathGraph])
%
% INPUT:
%   BSupTime          - time of the BSup pulses (180 degree inversion) in ms and before the start of readout, vector NBx1 (REQUIRED)
%   PresaturationTime - time in ms before the start of the readout, scalar, when the slice has been saturated (90 degree flip)
%                       this has to come before all the bSup pulses, but doesn't need to be always specified (OPTIONAL, DEFAULT = 0)
%   T1Time            - T1 time in ms of the observed tissue (OPTIONAL, DEFAULT = 1240 - GM at 3T)
%   SliceTime         - vector a time in ms when the slices are acquired after the start of the readout - typically, the first should
%                       be zero. Can be a vector of any order (OPTIONAL, DEFAULT = 0)
%   PathGraph         - If not empty, and a path is provided, then save a graph of the tissue saturation efficiency from presat to 
%                       the end of imaging. The extension should be JPG or will be replaced by JPG (OPTIONAL, DEFAULT = empty)
% OUTPUT:
%   SignalPercentage - signal percentage (0-1) of the remaining signal after background suppresion.
%                      it can be a scalar or a vector of the same dimension as sliceTime
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function assumes that the signal is, at first, optionally saturated by a 90 degree flip at PresaturationTime before readout.
%              Then follows a series of BSup pulses (times before readout are given) that do a 180 degree flip. The observed tissue relaxes with time
%              T1time and the signal attenuation is calculated for several slices acquired at times relative to the readout. 
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: signalPercentage = xASL_quant_BSupCalculation([500 100], 2700, 1240, [0 30 50 70],'/home/test/graph.jpg')
% __________________________________
% Copyright 2015-2020 ExploreASL

% Admin
if nargin < 1 || isempty(BSupTime)
	error('BSupTime is a required parameter');
end

if nargin < 2 || isempty(PresaturationTime)
	PresaturationTime = 0;
end

if numel(PresaturationTime) > 1
	error('PresaturationTime has to be a scalar');
end

if PresaturationTime
	if min(PresaturationTime > BSupTime) == 0
		error('PresaturationTime has to be longer than all BSupTimes');
	end
end

if nargin < 3 || isempty(T1Time)
	T1Time = 1240;
end

if nargin < 4 || isempty(SliceTime)
	SliceTime = 0;
end

if nargin < 5 || isempty(PathGraph)
	PathGraph = [];
end

if ~isempty(PathGraph)
	[fpath,fname,~] = fileparts(PathGraph);
	PathGraph = fullfile(fpath, [fname, '.jpg']);
end

% Sort the BSup pulses in descending order. Since in BIDS format the BSup pulses are given as time
% before the start of readout, descending order is from the first to last pulse
BSupTime = sort(BSupTime,'descend');
SignalPercentageInitial = zeros(1,1+numel(BSupTime));

% First phase is from the initial saturation to the first BSupPulse
if PresaturationTime
	% Recovery from 0 to 1 at T1 for the period between PresaturationTime and the first BSup pulse
	SignalPercentageInitial(1) = 1 - (1 - 0)*exp(-(PresaturationTime-BSupTime(1))/T1Time);
	PresaturationTimeDummy = PresaturationTime;
else
	% Without saturation, full signal is preset
	SignalPercentageInitial(1) = 1;
	% For the purpose of visualization, in the absence of presaturation we set the presaturation time to the 
	% maximal BSup
	PresaturationTimeDummy = max(BSupTime)+2;
end

% First phase calculation for a graph visual output
if ~isempty(PathGraph)
	% Preallocate
	SignalPercentageVector = zeros(PresaturationTimeDummy + max(SliceTime),1);
	if PresaturationTime
		SignalPercentageVector(1:PresaturationTimeDummy-BSupTime(1)-1,1) = 1 - (1 - 0)*exp(-(1:PresaturationTimeDummy-BSupTime(1)-1)/T1Time);
	else
		SignalPercentageVector(1:PresaturationTimeDummy-BSupTime(1)-1,1) = 1;
	end
end

% Iterate over the BSup pulses
for iBSup = 1:length(BSupTime)
	% First invert the signal by the BSuppulse
	SignalPercentageInitial(iBSup) = -SignalPercentageInitial(iBSup);
	
	% Then it relaxes
	if iBSup == length(BSupTime)
		% either until the readout
		SignalPercentageInitial(iBSup+1) = 1 - (1 - SignalPercentageInitial(iBSup))*exp(-BSupTime(iBSup)/T1Time);
	else
		% or until the next BSup pulse
		SignalPercentageInitial(iBSup+1) = 1 - (1 - SignalPercentageInitial(iBSup))*exp(-(BSupTime(iBSup)-BSupTime(iBSup+1))/T1Time);
	end
	
	% Repeat the same in a 1ms interval for a graphical output
	if ~isempty(PathGraph)
		% First invert the signal by the BSuppulse
		SignalPercentageVector(PresaturationTimeDummy-BSupTime(iBSup)) = -SignalPercentageVector(PresaturationTimeDummy-BSupTime(iBSup)-1);
		
		% Then it relaxes
		if iBSup == length(BSupTime)
			% either until the readout
			SignalPercentageVector(PresaturationTimeDummy-BSupTime(iBSup)+1:PresaturationTimeDummy) =...
				1 - (1 - SignalPercentageVector(PresaturationTimeDummy-BSupTime(iBSup)))*...
				exp(-(1:BSupTime(iBSup))/T1Time);
		else
			% or until the next BSup pulse
			SignalPercentageVector(PresaturationTimeDummy-BSupTime(iBSup)+1:PresaturationTimeDummy-BSupTime(iBSup+1)-1) =...
				1 - (1 - SignalPercentageVector(PresaturationTimeDummy-BSupTime(iBSup)))*...
				exp(-(1:BSupTime(iBSup)-BSupTime(iBSup+1)-1)/T1Time);
		end
	end
end

% SignalPercentageInitial is the signal percentage at the start of the readout

% Calculate the final signal percentages at the given slice times
SignalPercentage = 1 - (1 - SignalPercentageInitial(end))*exp(-(SliceTime)/T1Time);

% Repeat the same in a 1ms interval for a graphical output
if ~isempty(PathGraph)
	SignalPercentageVector(PresaturationTimeDummy+1:PresaturationTimeDummy+max(SliceTime)) = ...
		1 - (1 - SignalPercentageVector(PresaturationTimeDummy))*...
		exp(-(1:max(SliceTime))/T1Time);
end

% Plot the visualization
if ~isempty(PathGraph)
	% Create a figure without displaying it
	figureHandle = figure('visible','off');
	
	figureAxes = axes('Parent',figureHandle);
	
	% Plot the tissue signal percentage
	plot(SignalPercentageVector(:,1)*100);hold on
	
	% Plot a line at the start of the readout
	plot([PresaturationTimeDummy,PresaturationTimeDummy],[-100,100],'k')
	
	% Plot a vertical line at zero
	plot([1,PresaturationTimeDummy+max(SliceTime)],[0,0],'g');
	
	% Plot a horizontal line at each slice
	for iSliceTime = 2:length(SliceTime)
		plot(PresaturationTimeDummy+[SliceTime(iSliceTime),SliceTime(iSliceTime)],[-100,100],'r')
	end
	
	% Label the axes and time axis
	ylabel('Tissue signal percentage (%)');
	xlabel('Timing relative to the start of the readout (ms)');
	XTickLocation = [];
	XTickLabel    = {};
	% Tick at the start of labeling
	if PresaturationTime
		XTickLocation(1) = 0;
		XTickLabel{1}    = -PresaturationTime;
	end
	% Tick at each BSup inversion
	for iBSup = 1:length(BSupTime)
		XTickLocation(end+1) = PresaturationTimeDummy-BSupTime(iBSup);
		XTickLabel{end+1}    = num2str(-BSupTime(iBSup));
	end
	
	
	if length(SliceTime) > 1
		% Tick at second and last slice
		XTickLocation(end+1) = PresaturationTimeDummy+SliceTime(2);
		XTickLabel{end+1}    = num2str(SliceTime(2));
		if length(SliceTime) > 2
			XTickLocation(end+1) = PresaturationTimeDummy+SliceTime(end);
			XTickLabel{end+1}    = num2str(SliceTime(end));
		end
	else
		% Or at zero for 3D sequences and 2D sequences with a single slice
		XTickLocation(end+1) = PresaturationTimeDummy;
		XTickLabel{end+1}    = '0';
	end
	
	set(figureAxes,'XTick',XTickLocation,'XTickLabel',XTickLabel);
	
	% Save and close the figure
	saveas(figureHandle,PathGraph,'jpg');
	close(figureHandle);
end

end