function [DataOut] = xASL_str2num(DataIn)
%xASL_str2num str2num wrapper, replacing 'n/a' with NaN (BIDS convention)
% and converting only strings to numbers. Also allows inputting cells
%
% FORMAT:       [DataOut] = xASL_str2num(DataIn)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  str2num wrapper, replacing 'n/a' with NaN (BIDS convention)
%               and converting only strings to numbers. Also allows inputting cells.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL

if isempty(DataIn)
	DataOut = NaN;
elseif ~iscell(DataIn)
	DataOut = doConversion(DataIn);
else
	% For cells - do it one by one
	DataOut = zeros(size(DataIn));
	for ii = 1:length(DataIn)
		TempOut = doConversion(DataIn{ii});
		if isnumeric(TempOut)
			DataOut(ii) = TempOut;
		else
			DataOut(ii) = NaN;
		end
	end
end

return

function DataOut = doConversion(DataIn)
if isempty(DataIn)
    DataOut = NaN;
elseif ~isnumeric(DataIn)
        try % if non-numeric data, that isn't NaN,
            % then simply keep the data as is (if we cannot convert it)
            if strcmp(DataIn,'n/a') || strcmp(DataIn,'NaN')
                DataOut = NaN;
            else
                DataOut = str2num(DataIn);
            end
        catch
            DataOut = DataIn;
        end
else
    DataOut = DataIn;
end

if isempty(DataOut)
    DataOut = NaN;
end

return