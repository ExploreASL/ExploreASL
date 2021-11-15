function [resultText] = xASL_adm_BreakString(textToPrint, SymbolToFill, bColor, bNewLines, bPrintImmediately)
%xASL_adm_BreakString Pads symbols left and right of string and prints it
%
% FORMAT: [resultText] = xASL_adm_BreakString(textToPrint,SymbolToFill, bColor, bNewLines, bPrintImmediately)
% 
% INPUT:
%   textToPrint   - String containing a text (CHAR ARRAY, OPTIONAL, DEFAULT '')
%   SymbolToFill  - Symbol to fill the text (CHAR ARRAY, OPTIONAL, DEFAULT '=')
%   bColor        - Printf text in color if true (BOOLEAN, OPTIONAL, DEFAULT true)
%   bNewLines     - Add new lines
%                     - 0 - no new lines
%                     - 1 - a new line after 
%                     - 2 - new lines before and after (INTEGER, OPTIONAL, DEFAULT 2)
%   bPrintImmediately - Print the string immediately(BOOLEAN, OPTIONAL, DEFAULT true)
%
% OUTPUT:
%   resultText    - Padded string
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Pads symbols left and right of a string. By default it adds new lines, and colors, and prints the string, with
%                 a possibility to turn each of these options off.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        [resultText] = xASL_adm_BreakString('DCM2NII','=');
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% Define defaults
    if nargin<1 || isempty(textToPrint)
        textToPrint = '';
    end
    if nargin<2 || isempty(SymbolToFill)
        SymbolToFill = '=';
    end
    if nargin<3 || isempty(bColor)
        bColor = true;
    end
    if nargin<4 || isempty(bNewLines)
        bNewLines = 2;
    end
    if nargin<5
        bPrintImmediately = true;
    end

    % Define string width
    textWidth = 94;
    
    % Create the padded string
    resultText = xASL_adm_PaddedString(textToPrint, SymbolToFill, bColor, bNewLines, textWidth);
    
    % Print
    if bPrintImmediately
        fprintf(resultText);
    end


end


% Create the padded string
function resultText = xASL_adm_PaddedString(textToPrint, SymbolToFill, bColor, bNewLines, textWidth)

    %% Create default string
    resultText = repmat(SymbolToFill, 1, textWidth);
    
    % Pad spaces
    if ~isempty(textToPrint)
        textToPrint = [' ' textToPrint ' '];
    end

    % Check length
    if length(textToPrint)<100
        % Determine half of the length
        half = floor(length(textToPrint)/2);
        resultText(50-half:50-half+length(textToPrint)-1) = textToPrint;
    else
        fprintf(2,'Input string too long...\n');
    end
    
    %% Add specifics
    
    % Check if string should be colored
    if bColor
        resultText = sprintf(['[\b' resultText ']\b']);
    end
    
    % Add new lines
    if bNewLines==2
        resultText = sprintf(['\n' resultText '\n']);
	elseif bNewLines==1
        resultText = sprintf([resultText '\n']);
    end

end
