function [OutputNumber] = xASL_adm_CatchNumbersFromString(InputString)
%xASL_adm_CatchNumbersFromString Extracts a number from a char array.
%
% FORMAT: [OutputNumber] = xASL_adm_CatchNumbersFromString(InputString)
% 
% INPUT:
%   InputString   - String containing a number (CHAR ARRAY, REQUIRED)
%
% OUTPUT:
%   OutputNumber  - Number (DOUBLE)
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Extracts a number from a char array.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        [OutputNumber] = xASL_adm_CatchNumbersFromString('test123test');
% __________________________________
% Copyright 2015-2021 ExploreASL

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

