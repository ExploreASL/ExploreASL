function Nr = xASL_adm_ConvertTime2Nr(Time)
%xASL_adm_ConvertTime2Nr Converts time to number
% input hhmm -> output hh (with minutes in fractions/floating point)
% Inverse from xASL_adm_ConvertNr2Time
%
% FORMAT:       Nr = xASL_adm_ConvertTime2Nr(Time)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Converts time to number input hhmm -> output hh (with
%               minutes in fractions/floating point).
%               Inverse from xASL_adm_ConvertNr2Time.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________


if ~isnumeric( Time )
    Time    = str2num( Time );
end

HourN       = floor(Time./100);
MinN        = Time - (HourN.*100);

Nr          = HourN + (MinN./60);

end

