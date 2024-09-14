function stackText = xASL_test_StackToString(stack)
%xASL_test_StackToString Convert a stack of information to a string
%
% FORMAT: stackText = xASL_test_StackToString(stack)
%
% INPUT:        stack      - Stack (REQUIRED)
%
% OUTPUT:       stackText  - String
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Convert a stack of information to a string.
%
% EXAMPLE:      stackText = xASL_test_StackToString(stack);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

    % Fallback
    stackText = '';
    
    % Iterate over elements
    for iElement = 1:size(stack,1)
        thisFile = stack(iElement).file;
        thisName = stack(iElement).name;
        thisLine = stack(iElement).line;
        stackText = [stackText ', ' thisName ': line ' num2str(thisLine)];
    end
    
    % Remove initial ' ,'
    if ~isempty(stackText) && length(stackText)>3
        stackText = stackText(3:end);
    end
end
