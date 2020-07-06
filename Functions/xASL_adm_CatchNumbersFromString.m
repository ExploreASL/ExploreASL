function [OutputNumber] = xASL_adm_CatchNumbersFromString(InputString)
%xASL_adm_CatchNumbersFromString Summary of this function goes here
%   Detailed explanation goes here
%
% FORMAT: [OutputNumber] = xASL_adm_CatchNumbersFromString(InputString)
% 
% INPUT:
%   ...
%
% OUTPUT:
%   ...
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  ...
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL

    % Catch all numbers
    OutputNumber = NaN;
    B = regexp(InputString,'\d*','Match');
    for ii=1:length(B)
        if ~isempty(B{ii})
            OutputNumber(ii)=str2double(B{ii});
        else
            OutputNumber(ii)=NaN;
        end
    end 

end

