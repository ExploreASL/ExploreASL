function [DataOut] = xASL_num2str(DataIn, f,bConcatenate)
%xASL_num2str Wrapper around Matlab builtin 'num2str', bypassing strings/characters
% & BIDS-compatible
%
% FORMAT: [DataOut] = xASL_num2str(DataIn[, f,bConcatenate])
%
% INPUT:
%   DataIn 	- input data (can be any format) (REQUIRED)
%   f  - second argument of "num2str". Can be a format, but mostly used as number of characters to print (effectively rounding) (OPTIONAL)
%   bConcatenate - concatenate multiple-lines to a single line with ',' as a delimiter (OPTIONAL, DEFAULT = 0)
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
% xASL_num2str([1 15], 2,1)
% ans = '1,15'
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright (C) 2015-2019 ExploreASL
%
% 2019-05-02 HJM
% -----------------------------------------------------------------------------------------------------------------------------------------------------

if nargin < 2 || isempty(f)
	f = '';
end

if nargin < 3 || isempty(bConcatenate)
	bConcatenate = 0;
end

if isnumeric(DataIn)
	if isnan(DataIn)
		DataOut = 'n/a';
	elseif isempty(f)
		DataOut = num2str(DataIn);
	else
		DataOut = num2str(DataIn, f);
	end
	
	% If the bConcatenate option is ON, then check if concatenation is needed
	if bConcatenate
		if size(DataOut,1) > 1
			DataOut(1:end-1,end+1) = ',';
			DataOut = DataOut';
			DataOut = DataOut(:)';
		end
	end
else % bypass if not a number
    DataOut     = DataIn;
end


end

