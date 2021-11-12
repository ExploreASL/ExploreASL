function [resultText] = xASL_adm_PaddedString(textToPrint,SymbolToFill,Color,newLines,printImmediately)
%xASL_adm_PaddedString Pads symbols left and right of string
%
% FORMAT: [resultText] = xASL_adm_PaddedString(textToPrint,SymbolToFill)
% 
% INPUT:
%   textToPrint   - String containing a text (CHAR ARRAY, OPTIONAL)
%   SymbolToFill  - Symbol to fill the text (CHAR ARRAY, OPTIONAL)
%   Color         - Default or colored (BOOLEAN, OPTIONAL)
%   newLines      - Add new lines (BOOLEAN, OPTIONAL)
%   printImmediately - Print string (BOOLEAN, OPTIONAL)
%
% OUTPUT:
%   resultText    - Padded string
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Pads symbols left and right of string.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        [resultText] = xASL_adm_PaddedString('DCM2NII','=');
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% Define defaults
    if nargin<1
        textToPrint = '';
    end
    if nargin<2
        SymbolToFill = '=';
    end
    if nargin<3
        Color = true;
    end
    if nargin<4
        newLines = true;
    end
    if nargin<5
        printImmediately = true;
    end
    width = 94;
    
    %% Create default string
    resultText = repmat(SymbolToFill,1,width);
    
    % Pad spaces
    textToPrint = [' ' textToPrint ' '];
    
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
    if newLines
        resultText = sprintf(['\n' resultText '\n']);
    end
    
    % Print
    if printImmediately
        fprintf(resultText);
    end


end


