function SignalPercentage = xASL_quant_BSupCalculation(BSupTime, SaturationTime, T1Time, SliceTime)
% xASL_quant_BSupCalculation calculates the percentage of the remaining signal after a Bsup correction
%
% FORMAT: signalPercentage = xASL_quant_BSupCalculation(BSupTime[, SaturationTime, T1Time, SliceTime])
%
% INPUT:
%   BSupTime       - time of the BS pulses (180 degree inversion) in ms and before the start of readout, vector NBx1 (REQUIRED)
%   SaturationTime - time in ms before the start of the readout, scalar, when the slice has been saturated (90 degree flip)
%                    this has to come before all the bSup pulses, but doesn't need to be always specified (OPTIONAL, DEFAULT = 0)
%   T1Time         - T1 time in ms of the observed tissue (OPTIONAL, DEFAULT = 1240 - GM at 3T)
%   SliceTime      - vector a time in ms when the slices are acquired after the start of the readout - typically, the first should
%                    be zero. Can be a vector of any order (OPTIONAL, DEFAULT = 0)
% OUTPUT:
%   SignalPercentage - signal percentage (0-1) of the remaining signal after background suppresion.
%                      it can be a scalar or a vector of the same dimension as sliceTime
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function assumes that the signal is, at first, optionally saturated by a 90 degree flip at saturationTime before readout.
%              Then follows a series of BS pulses (times before readout are given) that do a 180 degree flip. The observed tissue relaxes with time
%              T1time and the signal attenuation is calculated for several slices acquired at times relative to the readout. 
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: signalPercentage = xASL_quant_BSupCalculation([500 100], 2700, 1240, [0 30 50 70])
% __________________________________
% Copyright 2015-2020 ExploreASL

% Admin
if nargin < 1 || isempty(BSupTime)
	error('BSupTime is a required parameter');
end

if nargin < 2 || isempty(SaturationTime)
	SaturationTime = 0;
end

if numel(SaturationTime) > 1
	error('SaturationTime has to be a scalar');
end

if SaturationTime
	if min(SaturationTime > BSupTime) == 0
		error('SaturationTime has to be longer than all BSupTimes');
	end
end

if nargin < 3 || isempty(T1Time)
	T1Time = 1240;
end

if nargin < 4 || isempty(SliceTime)
	SliceTime = 0;
end

% Sort the BS pulses descending - i.e. from the first to last
BSupTime = sort(BSupTime,'descend');
SignalPercentageInitial = zeros(1,1+numel(BSupTime));

% First phase is from the initial saturation to the first BSPulse
if SaturationTime
	% Recovery from 0 to 1 at T1 for the period between SaturationTime and the first BSpulse
	SignalPercentageInitial(1) = 1 - (1 - 0)*exp(-(SaturationTime-BSupTime(1))/T1Time);
else
	% Without saturation, full signal is preset
	SignalPercentageInitial(1) = 1;
end

% Iterate over the BSup pulses
for iBSup = 1:length(BSupTime)
	% First invert the signal by the BSpulse
	SignalPercentageInitial(iBSup) = -SignalPercentageInitial(iBSup);
	
	% Then it relaxes
	if iBSup == length(BSupTime)
		% either until the readout
		SignalPercentageInitial(iBSup+1) = 1 - (1 - SignalPercentageInitial(iBSup))*exp(-BSupTime(iBSup)/T1Time);
	else
		% or until the next BS pulse
		SignalPercentageInitial(iBSup+1) = 1 - (1 - SignalPercentageInitial(iBSup))*exp(-(BSupTime(iBSup)-BSupTime(iBSup+1))/T1Time);
	end
end

% SignalPercentageInitial is the signal percentage at the start of the readout

% Calculate the final signal percentages at the given slice times
SignalPercentage = 1 - (1 - SignalPercentageInitial(end))*exp(-(SliceTime)/T1Time);

end