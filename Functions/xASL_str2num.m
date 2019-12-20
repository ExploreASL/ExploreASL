function [DataOut] = xASL_str2num(DataIn)
%xASL_str2num str2num wrapper, replacing 'n/a' with NaN (BIDS convention)
% and converting only strings to numbers. Also allows inputting cells

if isempty(DataIn)
    DataOut = NaN;
elseif iscell(DataIn) && length(DataIn)>1
    error('Too long input');
elseif iscell(DataIn)
    DataIn = DataIn{1};
end

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

end