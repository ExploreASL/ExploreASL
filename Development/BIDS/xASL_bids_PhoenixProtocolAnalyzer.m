function [xasl,parameters] = xASL_bids_PhoenixProtocolAnalyzer(parameterList)
%xASL_bids_PhoenixProtocolAnalyzer This function analyzes the parameter list of the phoenix protocol (tag = [0x29,0x1020]).
%
% FORMAT: [xasl,parameters] = xASL_bids_PhoenixProtocolAnalyzer(parameterList);
%
% INPUT:
%        parameterList      - list of parameters from the reduced phoenix protocol (can be generated using xASL_bids_PhoenixProtocolReader) (REQUIRED)
%
% OUTPUT:
%        xasl               - list of xasl parameters
%        parameters         - list of phoenix parameters you want to know
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      This function analyzes the parameter list of the phoenix protocol (tag = [0x29,0x1020]).
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          [xasl,parameters] = xASL_bids_PhoenixProtocolAnalyzer(parameterList);
%
% REFERENCES:       ...
% __________________________________
% Copyright @ 2015-2020 ExploreASL


    %% Defaults
    parameters.tSequenceFileName = NaN;
    parameters.tProtocolName = NaN;
    
    
    %% Get the predefined parameters
    parameters = getIDofParameter(parameters,parameterList,false);
    
    
    %% Get xASL parameters
    
    if contains(parameters.tSequenceFileName,'SiemensSeq')
        xasl.SequenceType = 'Siemens'; 
    elseif contains(parameters.tSequenceFileName,'xxxxxxxxxxxxxx')
        xasl.SequenceType = 'Customer'; 
    end
    
    if contains(parameters.tSequenceFileName,'ep2d')
        xasl.Sequence = '2D EPI'; 
    elseif contains(parameters.tSequenceFileName,'xxxxxxxxxxxxxx')
        xasl.Sequence = '3D GRASE'; 
    end
    
    if contains(parameters.tSequenceFileName,'pasl')
       xasl.LabelingType = 'PASL'; 
    end
    


end


%% Get ID of parameter name
function parameters = getIDofParameter(parameters,phoenixParameterList,debugMode)

    % Get number of parameters
    parameterNames = fieldnames(parameters);
    numberOfParameters = numel(parameterNames);
    IDs = nan(numberOfParameters,1);
    values = nan(numberOfParameters,1);

    for curParameter=1:numberOfParameters
        % Get parameter name and ID in parameter list
        curName = parameterNames{curParameter,1};
        % Print out current information
        if debugMode
           fprintf('%s\n',curName); 
        end
        
        % Find parameter
        IDs(curParameter,1) = find(contains(phoenixParameterList, curName));
        if ~isempty(IDs(curParameter,1))
            parameters.(curName) = [phoenixParameterList{IDs(curParameter,1),2}];
        end
        
        
    end

    
    
end




