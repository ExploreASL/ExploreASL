function [stringValue] = xASL_adm_formatWithRoundedMagnitude(value, sigFigs)
% xASL_adm_formatWithRoundedMagnitude Formats a number with specified significant figures
%   
% FORMAT: str = xASL_adm_formatWithRoundedMagnitude(value, sigFigs)
%
% INPUT:
%   value               - numerical value to be formatted (REQUIRED)
%   sigFigs             - number of significant figures to format to (REQUIRED)
%
% OUTPUT: 
%   stringValue         - formatted string representing the number rounded to the specified significant figures
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  formatWithRoundedMagnitude takes a numerical input and formats it as a string with a fixed number of significant figures. 
%               The formatting avoids scientific notation, ensuring readability, while dynamically adjusting the number of decimal places 
%               based on the value's magnitude.
%
%               The function performs the following steps:
%               - Determines the order of magnitude of the input value.
%               - Calculates the scale factor to round the value to the desired significant figures.
%               - Dynamically adjusts the number of decimal places to maintain the specified significant figures.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:
%   stringValue = xASL_adm_formatWithRoundedMagnitude(0.05678000, 3);
%   % Result: '0.0568'
%   stringValue = xASL_adm_formatWithRoundedMagnitude(123400.0, 3);
%   % Result: '123000'
%
% __________________________________
% SPDX-License-Identifier: Apache-2.0
% __________________________________

    if value == 0
        stringValue = '0'; % Special case for zero
    else
        % Determine the order of magnitude
        orderOfMagnitude = floor(log10(abs(value)));

        % Scale the number to round to the desired significant figures
        scaleFactor = 10^(sigFigs - orderOfMagnitude - 1);
        roundedValue = round(value * scaleFactor) / scaleFactor;

        % Determine the number of decimal places based on the order of magnitude
        decimalPlaces = max(sigFigs - orderOfMagnitude - 1, 0);

        % Create a format string dynamically
        formatSpec = ['%.' num2str(decimalPlaces) 'f'];

        % Format the rounded number
        stringValue = sprintf(formatSpec, roundedValue);

        % Remove trailing floating point zeros
        if decimalPlaces>0
            stringValue = regexprep(stringValue, '\.?0+$', '');
        end
    end
end