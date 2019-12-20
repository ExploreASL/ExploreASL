function [OutputNumber] = xASL_adm_CatchNumbersFromString(InputString)
%xASL_adm_CatchNumbersFromString Summary of this function goes here
%   Detailed explanation goes here

    % Catch all numbers
    OutputNumber = NaN;
    B = regexp(InputString,'\d*','Match');
    for ii=1:length(B)
        if ~isempty(B{ii})
            OutputNumber(ii)=str2double(B{ii});
        else
            OutputNumber(ii)=NaN;
        end
    end 

end

