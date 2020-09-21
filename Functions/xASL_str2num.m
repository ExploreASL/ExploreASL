function [DataOut] = xASL_str2num(DataIn, bKeepCell, bReplaceNonNumerical)
%xASL_str2num str2num wrapper, replacing 'n/a' with NaN (BIDS convention)
% and converting only strings to numbers. Also allows inputting cells
%
% FORMAT: [DataOut] = xASL_str2num(DataIn[, bKeepCell, bReplaceNonNumerical])
% 
% INPUT: DataIn - character, numerical or cell array
%        bKeepCell - boolean specifying if cell array should be
%                    preserved & strings inside the cell converted to numbers (TRUE),
%                    or if the cell array should be converted to numerical (FALSE)
%                    Note that this parameter will only be used if DataIn
%                    was a cell array (otherwise forced to false)
%                    (OPTIONAL, DEFAULT = FALSE)
%        bReplaceNonNumerical - boolean specifying if a character or
%                               cell without a number should be converted to NaN (TRUE) or
%                               kept unchanged (FALSE). Note that this
%                               parameter will only be used if bKeepCell is true
%                               (OPTIONAL, DEFAULT = TRUE)
%               
%
% OUTPUT: DataOut   - converted numerical data
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: str2num wrapper, which only converts strings to numbers, and allows inputting cells.
%              Also, it replaces 'n/a' with NaN (BIDS convention). And it
%              has some other functionality as described in bKeepCell &
%              bReplaceNonNumerical above.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [5 6] = xASL_str2num('5 6')
%          [5 6 NaN] = xASL_str2num({'5' '6' 'Weird'})
%          {[5 6 NaN]} = xASL_str2num({'5' '6' 'Weird'}, 1)
%          {[5 6 'Weird']} = xASL_str2num({'5' '6' 'Weird'}, 1, 0)
% __________________________________
% Copyright 2015-2020 ExploreASL


%% -------------------------------------------------------
%% Admin
if nargin<2 || isempty(bKeepCell)
    bKeepCell = false;
elseif ~iscell(DataIn)
    bKeepCell = false; % only keep cell if DataIn was cell
end
if nargin<3 || isempty(bReplaceNonNumerical)
    bReplaceNonNumerical = true;
elseif ~bKeepCell
    bReplaceNonNumerical = true; % only keep original input for a cell
end

if isempty(DataIn)
	DataOut = NaN;
elseif ~iscell(DataIn)
	DataOut = doConversion(DataIn, bKeepCell, bReplaceNonNumerical);
else
    DataOut = cellfun(@(y) doConversion(y, bKeepCell, bReplaceNonNumerical), DataIn);
end


end

    
%% ==================================================================================
function DataOut = doConversion(DataIn, bKeepCell, bReplaceNonNumerical)

    
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

if isempty(DataOut) || ~isnumeric(DataOut)
    if bReplaceNonNumerical
        DataOut = NaN;
    else
        DataOut = DataIn;
    end
end

if bKeepCell
    DataOut = {DataOut};
end


end