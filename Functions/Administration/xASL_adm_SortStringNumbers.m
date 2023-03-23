function [outputString] = xASL_adm_SortStringNumbers(inputString, numberPatternInString)
%xASL_adm_SortStringNumbers Sort strings based on numbers
%
% FORMAT: [outputString] = xASL_adm_SortStringNumbers(inputString[, numberPatternInString])
%
% INPUT:
%   inputString             - cell array with strings that need numerical sorting (REQUIRED)
%   numberPatternInString   - descriptor of string format, with numerical
%                             location (OPTIONAL, DEFAULT = 'ASL_%d')
%
% OUTPUT:
%   outputString            - same as inputString but sorted numerically
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function sorts strings based on numbers
%               Normal matlab sort would sort {'ASL_1', 'ASL_2', 'ASL_12'}
%               to {'ASL_1', 'ASL_12', 'ASL_2'},
%               this function sorts them
%               to {'ASL_1', 'ASL_2', 'ASL_12'}
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: outputString = xASL_adm_SortStringNumbers(inputString);
% __________________________________
% Copyright 2015-2023 ExploreASL


    if nargin<2 || isempty(numberPatternInString)
        numberPatternInString = 'ASL_%d';
    elseif ~contains(numberPatternInString, '%d')
        error('numberPatternInString should contain a numerical descriptor %d');
    end

    % First sort alphabetically
    % (although numberPatternInString must support this)
    outputString = sort(inputString);
    
    % Then obtain numbers from text
    numberSorting = cellfun(@(x)sscanf(x, numberPatternInString), outputString);

    % sort them, and get the sorting order
      [~, sortOrder] = sort(numberSorting);
    % use to this sorting order to sort the filenames
      outputString = outputString(sortOrder);
end
