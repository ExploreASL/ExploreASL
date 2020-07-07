function [Nr DayInYear] = xASL_adm_ConvertDate2Nr(TempDate)
%xASL_adm_ConvertDate2Nr Converts date to number
% input mmdd -> output mm (with days in fractions/floating point)
% Inverse from ConvertNrDate
%
% FORMAT:       [Nr DayInYear] = xASL_adm_ConvertDate2Nr(TempDate)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Converts date to number input mmdd -> output mm (with days in fractions/floating point).
%               Inverse from ConvertNrDate.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL

DaysInMonth     = [31 28.25 31 30 31 30 31 31 30 31 30 31];

for iI=1:size(TempDate,1)
    for iJ=1:size(TempDate,2)
        for iK=1:size(TempDate,3)
            clear Date MonthN DayN
            
            Date    = TempDate(iI,iJ,iK);

            if ~isnumeric( Date )
                Date    = str2num( Date );
            end

            if  length(num2str(Date))==8 % remove the year
                Date    = num2str(Date);
                Date    = Date(5:8);
                Date    = str2num( Date );
            end

            MonthN      = floor(Date/100);
            DayN        = Date - (MonthN*100);

            DaysMonthCount  = sum(DaysInMonth(1:MonthN));
            DayInYear(iI,iJ,iK)       = DaysMonthCount + DayN;
            
            Nr(iI,iJ,iK)          = MonthN + (DayN/DaysInMonth(MonthN));
        end
    end
end

            

end

