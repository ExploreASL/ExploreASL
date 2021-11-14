function [resultText] = xASL_adm_BreakString(textToPrint,SymbolToFill,Color,newLines,printImmediately)
%xASL_adm_BreakString Pads symbols left and right of string and prints it
%
% FORMAT: [resultText] = xASL_adm_BreakString(textToPrint,SymbolToFill)
% 
% INPUT:
%   textToPrint   - String containing a text (CHAR ARRAY, OPTIONAL)
%   SymbolToFill  - Symbol to fill the text (CHAR ARRAY, OPTIONAL)
%   Color         - Default or colored (BOOLEAN, OPTIONAL)
%   newLines      - Add new lines, where 1 only means after and 2 means before and after (INTEGER, OPTIONAL)
%   printImmediately - Print string (BOOLEAN, OPTIONAL)
%
% OUTPUT:
%   resultText    - Padded string
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Pads symbols left and right of string.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        [resultText] = xASL_adm_BreakString('DCM2NII','=');
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% Define defaults
    if nargin<1
        textToPrint = '';
    end
    if nargin<2 || isempty(SymbolToFill)
        SymbolToFill = '=';
    end
    if nargin<3 || isempty(Color)
        Color = true;
    end
    if nargin<4 || isempty(newLines)
        newLines = 2;
    end
    if nargin<5
        printImmediately = true;
    end

    % Define string width
    textWidth = 94;
    
    % Create the padded string
    resultText = xASL_adm_PaddedString(textToPrint,SymbolToFill,Color,newLines,textWidth);
    
    % Print
    if printImmediately
        fprintf(resultText);
    end


end


% Create the padded string
function resultText = xASL_adm_PaddedString(textToPrint, SymbolToFill, Color, newLines, textWidth)

    %% Create default string
    resultText = repmat(SymbolToFill,1,textWidth);
    
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
    if Color
        resultText = sprintf(['[\b' resultText ']\b']);
    end
    
    % Add new lines
    if newLines==2
        resultText = sprintf(['\n' resultText '\n']);
    else
        resultText = sprintf([resultText '\n']);
    end


end