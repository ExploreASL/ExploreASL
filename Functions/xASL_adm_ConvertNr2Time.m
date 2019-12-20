function Time = xASL_adm_ConvertNr2Time( Nr )
%xASL_adm_ConvertNr2Time Converts number to time
% input hh (with minutes in fractions/floating point) -> output hhmm 
% Inverse from xASL_adm_ConvertTime2Nr

if ~isnumeric( Nr )
    Nr    = str2num( Nr );
end

HourN       = floor(Nr);
MinN        = round((Nr - HourN) * 60);

Time        = HourN*100 + MinN;

% add leading zeros
Time        = sprintf( '%.4d', Time );

end

