function TimeString = xASL_adm_ConvertSeconds2TimeString(Seconds)
%xASL_adm_ConvertNr2Time Converts number to time
% input hh (with minutes in fractions/floating point) -> output hhmm 
% Inverse from xASL_adm_ConvertTime2Nr
%
% FORMAT: TimeString = xASL_adm_ConvertSeconds2TimeString(Seconds)
%
% INPUT:
%   Seconds	    - Number
%
% OUTPUT:
%   TimeString  - Time
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Converts number to time
% input hh (with minutes in fractions/floating point) -> output hhmm 
% Inverse from xASL_adm_ConvertTime2Nr.
%
% EXAMPLE:     n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright 2015-2021 ExploreASL

if ~isnumeric(Seconds)
    error('Input parameter Seconds should be numeric');
end

HourN = floor(Seconds/3600);
Seconds = Seconds - HourN*3600;
MinN = floor(Seconds/60);
Seconds = Seconds - MinN*60;

if HourN>0
    TimeString = [sprintf('%.1d', HourN) 'h' sprintf('%.2d', MinN) 'm' sprintf('%.2d', Seconds)];
elseif MinN>0
    TimeString = [sprintf('%.2d', MinN) 'm' sprintf('%.2d', Seconds)];
else
    TimeString = [sprintf('%.2d', Seconds)];
end

end

