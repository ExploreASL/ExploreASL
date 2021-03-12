function [DataOut] = xASL_num2str(DataIn, f, bConcatenate, strDelimiter)
%xASL_num2str Wrapper around Matlab builtin 'num2str', bypassing strings/characters
% & BIDS-compatible
%
% FORMAT: [DataOut] = xASL_num2str(DataIn[, f, bConcatenate, strDelimiter])
%
% INPUT:
%   DataIn 	     - input data (can be any format) (REQUIRED)
%   f            - second argument of "num2str". Can be a format, but mostly used as number of characters to print (effectively rounding) (OPTIONAL)
%   bConcatenate - concatenate multiple-lines to a single line with ',' as a delimiter (OPTIONAL, DEFAULT = 0)
%   strDelimiter - string to be used as a delimiter for concatenating when bConcatenate is TRUE, empty delimiters are allowed (OPTIONAL, DEFAULT = ',')
%
% OUTPUT:
%   DataOut - numbers converted to string, or bypassed data
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: When the provided data is numeric, this function will convert the number to string/characters,
% rewriting NaN into 'n/a' (BIDS convention) but otherwise preserving the Matlab builtin functionality, also for the second argument "f".
% If non-numeric data is provided, it is bypassed (avoiding any issues "num2str" will have with non-numeric data).
% It can concatenate an array/matrix of strings, taking first the columns in the first row, and then going across the rows.
% See builtin num2str for more details
%   
% EXAMPLE: 
% xASL_num2str(10.5798)
% ans = '10.5798'
% xASL_num2str(10.5798, 2)
% ans = 11
% xASL_num2str([1;15], 2,1)
% ans = '1,15 '
% xASL_num2str([1,15;2,13], 2,1,',,')
% ans = '1,,15,,2,,13  '
% [DataOut] = xASL_num2str(123.456789, '%.5f'); % Get exactly 5 digits after the comma
% DataOut = '123.45679'
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright (C) 2015-2020 ExploreASL
%
% 2019-05-02 HJM
% -----------------------------------------------------------------------------------------------------------------------------------------------------

% Default: no format
if nargin < 2 || isempty(f)
	f = '';
end

% Default: set concatenate to 0
if nargin < 3 || isempty(bConcatenate)
	bConcatenate = 0;
end

% Default: set delimiter
if nargin < 4
	strDelimiter = ',';
elseif isempty(strDelimiter)
	strDelimiter = '';
end

if size(strDelimiter,1) > 1
	error('The delimiter can only have a single row');
end

if isnumeric(DataIn)
	if isnan(DataIn)
		DataOut = 'n/a';
	elseif isempty(f)
		DataOut = num2str(DataIn);
    elseif strcmp(f,'auto') % Automatic mode
        if floor(DataIn)==DataIn % Integers
            DataOut = num2str(DataIn, '%d');
        else
            DataOut = num2str(DataIn, '%.12f');
            % Remove trailing zeros
            DataOut = strip(DataOut,'right','0');
        end
	else
		DataOut = num2str(DataIn, f);
	end
	
	% If the bConcatenate option is ON, then check if concatenation is needed
	if bConcatenate
		if size(DataOut,1) > 1
			DataOut(1:end-1,(end+1):(end+length(strDelimiter))) = repmat(strDelimiter,[(size(DataOut,1)-1) 1]);
			DataOut = DataOut';
			DataOut = DataOut(:)';
		end
		% Replace the extra spaces with the delimiter
		[indStart,indEnd] = regexp(DataOut,'\d{1}\ +');
		for iSpace = 1:length(indStart)
			DataOut = [DataOut(1:indStart(iSpace)) strDelimiter DataOut(indEnd(iSpace)+1:end)];
			indShift = length(strDelimiter) - (indEnd(iSpace)-indStart(iSpace));
			indStart = indStart + indShift;
			indEnd   = indEnd + indShift;
		end
	end
else % bypass if not a number
    DataOut     = DataIn;
end


end

