function [DataOut] = xASL_num2str(DataIn, stringFormat, bConcatenate, stringDelimiter)
%xASL_num2str Wrapper around Matlab builtin 'num2str', bypassing strings/characters & BIDS-compatible
%
% FORMAT: [DataOut] = xASL_num2str(DataIn[, stringFormat, bConcatenate, stringDelimiter])
%
% INPUT:
%   DataIn 	        - Input data (can be any format) (REQUIRED, INTEGER or DOUBLE, SCALAR or ARRAY)
%   stringFormat    - Second argument of "num2str". Can be a format, but mostly used as number
%                     of characters to print (effectively rounding) (OPTIONAL, DEFAULT = '')
%   bConcatenate    - Concatenate multiple-lines to a single line with ',' as a delimiter (OPTIONAL, DEFAULT = 1)
%   stringDelimiter - String to be used as a delimiter for concatenating when bConcatenate 
%                     is TRUE, empty delimiters are allowed (OPTIONAL, DEFAULT = ',')
%
% OUTPUT:
%   DataOut         - Numbers converted to string, or bypassed data (CHAR ARRAY)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: When the provided data is numeric, this function will convert the number to 
%              string/characters, rewriting NaN into 'n/a' (BIDS convention) but otherwise 
%              preserving the Matlab builtin functionality, also for the second argument "f".
%              If non-numeric data is provided, it is bypassed (avoiding any issues "num2str"
%              will have with non-numeric data).
%              It can concatenate an array/matrix of strings, taking first the columns in the 
%              first row, and then going across the rows. See builtin num2str for more details. 
%              Column vectors are converted to row vectors.
%   
% EXAMPLE:
%
% xASL_num2str(10.5798)
% ans = '10.5798'
% xASL_num2str(10.5798, 2)
% ans = '11'
% xASL_num2str([1;15], 2,1)
% ans = '1,15'
% xASL_num2str([1,15;2,13], 2,1,',,')
% ans = '1,,15,,2,,13'
% Get exactly 5 digits after the comma
% DataOut = xASL_num2str(123.456789, '%.5f');
% DataOut = '123.45679'
%
% Automatic mode (remove trailing zeros)
% DataOut = xASL_num2str(1.23456789000, 'auto')
% DataOut = '1.23456789'
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright (c) 2015-2021 ExploreASL

% Check input
if isnumeric(DataIn)
    % Convert row to column vector on default (for later concatenation and delimiter)
    if size(DataIn,2)>size(DataIn,1)
        DataIn = DataIn';
    end
end

% Default: no format
if nargin < 2 || isempty(stringFormat)
	stringFormat = '';
end

% Default: set concatenate to 1
if nargin < 3 || isempty(bConcatenate)
	bConcatenate = 1;
end

% Default: set delimiter
if nargin < 4
	stringDelimiter = ',';
elseif isempty(stringDelimiter)
	stringDelimiter = '';
end

if size(stringDelimiter,1) > 1
	error('The delimiter can only have a single row');
end

% Do the conversion if DataIn is numeric
if isnumeric(DataIn)
	if isnan(DataIn)
		DataOut = 'n/a';
	elseif isempty(stringFormat)
		DataOut = num2str(DataIn);
    elseif strcmp(stringFormat,'auto') % Automatic mode
        if floor(DataIn)==DataIn % Integers
            DataOut = num2str(DataIn, '%d');
        else
            % Convert to specific format
            DataOut = num2str(DataIn, '%.12f');
            % Remove trailing zeros
            DataOut = xASL_adm_RemoveTrailingSymbol(DataOut,'0');
            % Fix perfect integers (previously 2.0000 was converted to 2. instead of 2.0)
            if strcmp(DataOut(end),'.')
                DataOut = [DataOut '0'];
            end
        end
	else
		DataOut = num2str(DataIn, stringFormat);
	end
	
    % If the bConcatenate option is ON, then check if concatenation is needed
    if bConcatenate
        if size(DataOut,1) > 1
            DataOut(1:end-1,(end+1):(end+length(stringDelimiter))) = repmat(stringDelimiter,[(size(DataOut,1)-1) 1]);
            DataOut = DataOut';
            DataOut = DataOut(:)';
        end
        % Replace the extra spaces with the delimiter
        [indStart,indEnd] = regexp(DataOut,'\d{1}\ +');
        for iSpace = 1:length(indStart)
            DataOut = [DataOut(1:indStart(iSpace)) stringDelimiter DataOut(indEnd(iSpace)+1:end)];
            indShift = length(stringDelimiter) - (indEnd(iSpace)-indStart(iSpace));
            indStart = indStart + indShift;
            indEnd   = indEnd + indShift;
        end
    end
    
    % Always convert back to row afterwards
    DataOut = DataOut(:)';
    
    % Remove trailing spaces
    DataOut = xASL_adm_RemoveTrailingSymbol(DataOut);
    
else % bypass if not a number
    DataOut = DataIn;
end


end


%% Strip trailing symbol from string
function stringToStrip = xASL_adm_RemoveTrailingSymbol(stringToStrip,removeMe)

    % Default: remove spaces
    if nargin<2
        removeMe = 'wspace';
    end

    % Remove trailing symbols in Matlab versions > 2016b
    [~, dateString] = version();
    if str2num(dateString(end-3:end))>2016 && ~strcmp(removeMe,'wspace')
        stringToStrip = strip(stringToStrip,'right',removeMe);
    elseif strcmp(removeMe,'wspace')
        indToStrip = isstrprop(stringToStrip,'wspace');
        stringToStrip(indToStrip) = [];
		% And also remove zero string
		indToStrip = stringToStrip == 0;
		stringToStrip(indToStrip) = [];
    else
        % Find out trailing zeros
        indSymbol = regexp(stringToStrip,[removeMe '*$']);
        % If the string ends with those symbols only
        if ~isempty(indSymbol)
            % Remove the trailing symbols
            stringToStrip = stringToStrip(1:(indSymbol(1)-1));
        end
    end

end



