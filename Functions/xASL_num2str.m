function [DataOut] = xASL_num2str(DataIn, f)
%xASL_num2str Wrapper around Matlab builtin 'num2str', bypassing strings/characters
% & BIDS-compatible
%
% FORMAT: [DataOut] = xASL_num2str(DataIn[, f])
%
% INPUT:
%   DataIn 	- input data (can be any format) (REQUIRED)
%   f  - second argument of "num2str". Can be a format, but mostly used as number of characters to print (effectively rounding) (OPTIONAL)
%
% OUTPUT:
%   DataOut - numbers converted to string, or bypassed data
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: when the provided data is numeric, this function will convert the number to string/characters,
% rewriting NaN into 'n/a' (BIDS convention) but otherwise preserving the Matlab builtin functionality, also for the second argument "f".
% If non-numeric data is provided, it is bypassed (avoiding any issues "num2str" will have with non-numeric data).
% See builtin num2str for more details
%   
% EXAMPLE: 
% xASL_num2str(10.5798)
% ans = '10.5798'
% xASL_num2str(10.5798, 2)
% ans = 11
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright (C) 2015-2019 ExploreASL
%
% 2019-05-02 HJM
% -----------------------------------------------------------------------------------------------------------------------------------------------------

if isnumeric(DataIn)
    if isnan(DataIn)
        DataOut = 'n/a';
    elseif nargin<2
        DataOut = num2str(DataIn);
    else
        DataOut = num2str(DataIn, f);
    end
else % bypass if not a number
    DataOut     = DataIn;
end


end

