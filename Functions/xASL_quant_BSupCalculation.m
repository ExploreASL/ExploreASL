function SignalPercentage = xASL_quant_BSupCalculation(BackgroundSuppressionPulseTime, PresaturationTime, T1Time, SliceTime, PathGraph, bFigureExists)
% xASL_quant_BSupCalculation calculates the percentage of the remaining signal after a Bsup correction
%
% FORMAT: signalPercentage = xASL_quant_BSupCalculation(BackgroundSuppressionPulseTime[, PresaturationTime, T1Time, SliceTime, PathGraph])
%
% INPUT:
%   BackgroundSuppressionPulseTime          - time of the BSup pulses (180 degree inversion) in ms and before the start of readout, vector NBx1 (REQUIRED)
%   PresaturationTime - time in ms before the start of the readout, scalar, when the slice has been saturated (90 degree flip)
%                       this has to come before all the bSup pulses, but doesn't need to be always specified (OPTIONAL, DEFAULT = 0)
%   T1Time            - T1 time in ms of the observed tissue (OPTIONAL, DEFAULT = 1240 - GM at 3T)
%   SliceTime         - vector a time in ms when the slices are acquired after the start of the readout - typically, the first should
%                       be zero. Can be a vector of any order (OPTIONAL, DEFAULT = 0)
%   PathGraph         - If not empty, and a path is provided, then save a graph of the tissue signal percentage from the presaturation pulse to 
%                       the end of imaging. The extension will be set to JPG (OPTIONAL, DEFAULT = empty)
%   bFigureExists - If a figure handle exists, the graph will be
%                        created but not saved to disk (e.g. for
%                        subplotting inside ExploreASL's xASL_quant_M0.m)
% OUTPUT:
%   SignalPercentage - signal percentage (0-1) of the remaining signal after background suppression.
%                      It can be a scalar or a vector of the same dimension as sliceTime
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function computes the tissue signal percentage that
%              remains after background suppression pulses are played in the ASL
%              acquisition.
%              It assumes that the signal is, at first, optionally saturated by a 90 degree flip at PresaturationTime before readout.
%              Then follows a series of BSup pulses (times before readout are given) that do a 180 degree flip. The observed tissue relaxes with time
%              T1time and the signal attenuation is calculated for several slices acquired at times relative to the readout. 
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: signalPercentage = xASL_quant_BSupCalculation([500 100], 2700, 1240, [0 30 50 70], '/home/test/graph.jpg')
% __________________________________
% Copyright 2015-2020 ExploreASL
 
% Admin
if nargin < 1 || isempty(BackgroundSuppressionPulseTime)
    error('BackgroundSuppressionPulseTime is a required parameter');
end
 
if nargin < 2 || isempty(PresaturationTime)
    PresaturationTime = 0;
end
 
if numel(PresaturationTime) > 1
    error('PresaturationTime has to be a scalar');
end
 
if PresaturationTime
    if min(PresaturationTime > BackgroundSuppressionPulseTime) == 0
        error('PresaturationTime has to be longer than all BackgroundSuppressionPulseTimes');
    end
end
 
% Check order of magnitude
if PresaturationTime<100
    warning(['PresaturationTime=' num2str(PresaturationTime) ', seems invalid']);
end
if BackgroundSuppressionPulseTime<5
    warning(['BackgroundSuppressionTime of ' num2str(BackgroundSuppressionPulseTime) ', seems invalid']);
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
 
if nargin < 6 || isempty(bFigureExists)
    bFigureExists = 0;
end

if ~isempty(PathGraph)
    [fpath, fname] = fileparts(PathGraph);
    PathGraph = fullfile(fpath, [fname '.jpg']);
end
 
% Sort the BSup pulses in descending order. Since in BIDS format the BSup pulses are given as time
% before the start of readout, descending order is from the first to last pulse
BackgroundSuppressionPulseTime = sort(BackgroundSuppressionPulseTime, 'descend');
SignalPercentageInitial = zeros(1, 1 + numel(BackgroundSuppressionPulseTime));
 
% First phase is from the initial saturation to the first BSupPulse
if PresaturationTime
    % Recovery from 0 to 1 at T1 for the period between PresaturationTime and the first BSup pulse
    SignalPercentageInitial(1) = 1 - (1 - 0)*exp(-(PresaturationTime-BackgroundSuppressionPulseTime(1))/T1Time);
    PresaturationTimeDummy = PresaturationTime;
else
    % Without saturation, full signal is preset
    SignalPercentageInitial(1) = 1;
    % For the purpose of visualization, in the absence of presaturation we set the presaturation time to the 
    % maximal BSup
    PresaturationTimeDummy = max(BackgroundSuppressionPulseTime)+2;
end
 
% First phase calculation for a graph visual output
if ~isempty(PathGraph)
    % Preallocate
    SignalPercentageVector = zeros(round(PresaturationTimeDummy + max(SliceTime)), 1);
    if PresaturationTime
        SignalPercentageVector(1:PresaturationTimeDummy-BackgroundSuppressionPulseTime(1)-1,1) = 1 - (1 - 0)*exp(-(1:PresaturationTimeDummy-BackgroundSuppressionPulseTime(1)-1)/T1Time);
    else
        SignalPercentageVector(1:PresaturationTimeDummy-BackgroundSuppressionPulseTime(1)-1,1) = 1;
    end
end
 
% Iterate over the BSup pulses
for iBSup = 1:length(BackgroundSuppressionPulseTime)
    % First invert the signal by the BSuppulse
    SignalPercentageInitial(iBSup) = -SignalPercentageInitial(iBSup);
    
    % Then it relaxes
    if iBSup == length(BackgroundSuppressionPulseTime)
        % either until the readout
        SignalPercentageInitial(iBSup+1) = 1 - (1 - SignalPercentageInitial(iBSup))*exp(-BackgroundSuppressionPulseTime(iBSup)/T1Time);
    else
        % or until the next BSup pulse
        SignalPercentageInitial(iBSup+1) = 1 - (1 - SignalPercentageInitial(iBSup))*exp(-(BackgroundSuppressionPulseTime(iBSup)-BackgroundSuppressionPulseTime(iBSup+1))/T1Time);
    end
    
    % Repeat the same in a 1ms interval for a graphical output
    if ~isempty(PathGraph)
        % First invert the signal by the BSuppulse
        SignalPercentageVector(PresaturationTimeDummy-BackgroundSuppressionPulseTime(iBSup)) = -SignalPercentageVector(PresaturationTimeDummy-BackgroundSuppressionPulseTime(iBSup)-1);
        
        % Then it relaxes
        if iBSup == length(BackgroundSuppressionPulseTime)
            % either until the readout
            SignalPercentageVector(PresaturationTimeDummy-BackgroundSuppressionPulseTime(iBSup)+1:PresaturationTimeDummy) =...
                1 - (1 - SignalPercentageVector(PresaturationTimeDummy-BackgroundSuppressionPulseTime(iBSup)))*...
                exp(-(1:BackgroundSuppressionPulseTime(iBSup))/T1Time);
        else
            % or until the next BSup pulse
            SignalPercentageVector(PresaturationTimeDummy-BackgroundSuppressionPulseTime(iBSup)+1:PresaturationTimeDummy-BackgroundSuppressionPulseTime(iBSup+1)-1) =...
                1 - (1 - SignalPercentageVector(PresaturationTimeDummy-BackgroundSuppressionPulseTime(iBSup)))*...
                exp(-(1:BackgroundSuppressionPulseTime(iBSup)-BackgroundSuppressionPulseTime(iBSup+1)-1)/T1Time);
        end
    end
end
 
% SignalPercentageInitial is the signal percentage at the start of the readout
 
% Calculate the final signal percentages at the given slice times
SignalPercentage = 1 - (1 - SignalPercentageInitial(end))*exp(-(SliceTime)/T1Time);
 
% Repeat the same in a 1ms interval for a graphical output
if ~isempty(PathGraph)
    SignalPercentageVector(PresaturationTimeDummy+1:round(PresaturationTimeDummy+max(SliceTime))) = ...
        1 - (1 - SignalPercentageVector(PresaturationTimeDummy))*...
        exp(-(1:max(SliceTime))/T1Time);
end
 
% Plot the visualization
if ~isempty(PathGraph)
    % Create a figure without displaying it
    if ~bFigureExists
        figureHandle = figure('visible','off');
        figureAxes = axes('Parent', figureHandle);
    end
    
    % Plot the tissue signal percentage
    plot(SignalPercentageVector(:,1)*100); hold on
    
    % Plot a line at the start of the readout
    plot([PresaturationTimeDummy,PresaturationTimeDummy],[-100,100],'g')
    
    % Plot a vertical line at zero
    plot([1,PresaturationTimeDummy+max(SliceTime)],[0,0],'k');
    
    % Plot a horizontal line at each slice
    for iSliceTime = 2:length(SliceTime)
        plot(PresaturationTimeDummy+[SliceTime(iSliceTime),SliceTime(iSliceTime)],[-100,100],'r')
    end
    
    % Label the axes and time axis
    title('Background suppression effect on tissue signal');
    ylabel('Tissue signal percentage (%)');
    xlabel('Timing relative to the start of the readout (ms)');
    XTickLocation = [];
    XTickLabel    = {};
    % Tick at the start of labeling
    if PresaturationTime
        XTickLocation(1) = 0;
        XTickLabel{1} = -PresaturationTime;
    end
    % Tick at each BSup inversion
    for iBSup = 1:length(BackgroundSuppressionPulseTime)
        XTickLocation(end+1) = PresaturationTimeDummy-BackgroundSuppressionPulseTime(iBSup);
        XTickLabel{end+1} = num2str(-BackgroundSuppressionPulseTime(iBSup));
    end
    
    if length(SliceTime) > 1
        % Tick at second and last slice
        XTickLocation(end+1) = PresaturationTimeDummy+SliceTime(2);
        XTickLabel{end+1}    = num2str(SliceTime(2));
        if length(SliceTime) > 2
            XTickLocation(end+1) = PresaturationTimeDummy+SliceTime(end);
            XTickLabel{end+1} = num2str(SliceTime(end));
        end
    else
        % Or at zero for 3D sequences and 2D sequences with a single slice
        XTickLocation(end+1) = PresaturationTimeDummy;
        XTickLabel{end+1}    = '0';
    end
    
    if ~bFigureExists
        set(figureAxes, 'XTick', XTickLocation, 'XTickLabel', XTickLabel);
    
        % Save and close the figure
        saveas(figureHandle, PathGraph, 'jpg');
        close(figureHandle);
    end
end
 

end