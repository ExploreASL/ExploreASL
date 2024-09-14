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
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

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
