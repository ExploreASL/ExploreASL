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
    parameters{1,1} = 'tSequenceFileName';                      parameters{1,2} = NaN;
    parameters{2,1} = 'tProtocolName';                          parameters{2,2} = NaN;
    parameters{3,1} = 'UserScaleFactor';                        parameters{3,2} = NaN;
    parameters{4,1} = 'M0Threshold';                            parameters{4,2} = NaN;
    parameters{5,1} = 'TI1_us';                                 parameters{5,2} = NaN;
    parameters{6,1} = 'TI2_us';                                 parameters{6,2} = NaN;
    parameters{7,1} = 'sProtConsistencyInfo.tBaselineString';   parameters{7,2} = NaN;
    parameters{8,1} = 'sAsl.ulMode';                            parameters{8,2} = NaN;
    parameters{9,1} = 'alTR[0]';                                parameters{9,2} = NaN;
    parameters{10,1} = 'alTI[0]';                               parameters{10,2} = NaN;
    parameters{11,1} = 'alTI[1]';                               parameters{11,2} = NaN;
    parameters{12,1} = 'alTI[2]';                               parameters{12,2} = NaN;
    parameters{13,1} = 'alTE[0]';                               parameters{13,2} = NaN;
    parameters{14,1} = 'acFlowComp[0]';                         parameters{14,2} = NaN;
    parameters{15,1} = 'sGroupArray.sPSat.dThickness';          parameters{15,2} = NaN;
    parameters{16,1} = 'sGroupArray.sPSat.dGap';                parameters{16,2} = NaN;
    parameters{17,1} = 'lRepetitions';                          parameters{17,2} = NaN;
    parameters{18,1} = 'sAsl.fFlowLimit';                       parameters{18,2} = NaN;
    parameters{19,1} = 'sWipMemBlock.alFree[0]';                parameters{19,2} = NaN;
    
    %% Get the predefined parameters
    parameters = getPhoenixParameters(parameters,parameterList,false);
    
    %% Get xASL parameters
    parameters = convertCellArrayToStruct(parameters);

    %% Analyze parameters
    if contains(parameters.tSequenceFileName,'SiemensSeq')
        xasl.SequenceType = 'Siemens'; 
    elseif contains(parameters.tSequenceFileName,'CustomerSeq')
        xasl.SequenceType = 'Customer'; 
    end
    
    if contains(lower(parameters.tSequenceFileName),'ep2d')
        xasl.Sequence = '2D EPI'; 
    elseif contains(lower(parameters.tSequenceFileName),'gse')
        xasl.Sequence = '3D GRASE'; 
    end
    
    if contains(lower(parameters.tSequenceFileName),'pasl')
        xasl.LabelingType = 'PASL';
    end
    if contains(lower(parameters.tSequenceFileName),'casl')
        xasl.LabelingType = 'CASL';
    end
    if contains(lower(parameters.tSequenceFileName),'pcasl')
        xasl.LabelingType = 'PCASL';
    end
    
    if parameters.sAslulMode==4
        xasl.Technique = 'Q2TIPS';
    end
    
    if parameters.alTI0~=100000 && parameters.alTI2~=7000000
        xasl.ScanType = 'ASL';
    elseif parameters.alTI0==100000 && parameters.alTI2==7000000
        xasl.ScanType = 'PseudoM0';
    end



end

%% Get individual phoenix parameter
function value = getPhoePar(parameters,curParToExtract)
    % Get cell ID
    cellID = find(contains(parameters, curParToExtract));
    % Return value
    value = parameters{cellID,2};
    % Convert string numbers to numbers
    if ~isnan(str2double(value))
        value = str2double(value);
    end
end

%% Convert cell array to struct
function parameterStruct = convertCellArrayToStruct(parameters)
    % Create struct
    parameterStruct = struct;
    % Get number of parameters
    parameterNames = parameters(:,1);
    numberOfParameters = numel(parameterNames);
    % Convert each line to a field
    for curParameter=1:numberOfParameters
        curParToExtract = parameters{curParameter,1};
        value = getPhoePar(parameters,curParToExtract);
        validFieldName = matlab.lang.makeValidName(curParToExtract,'ReplacementStyle','delete');
        parameterStruct.(validFieldName) = value;        
    end
end

%% Get ID of parameter name
function parameters = getPhoenixParameters(parameters,phoenixParameterList,debugMode)

    % Get number of parameters
    parameterNames = parameters(:,1);
    numberOfParameters = numel(parameterNames);
    IDs = nan(numberOfParameters,1);
    newPar = 1;

    for curParameter=1:numberOfParameters
        % Get parameter name and ID in parameter list
        curName = parameterNames{curParameter,1};
        
        % Find parameter
        curFoundParID = find(contains(phoenixParameterList, curName));
        if length(curFoundParID)==1
            IDs(curParameter,1) = curFoundParID;
        elseif length(curFoundParID)>1            
            % Select exact match
            for tmpPar=1:length(curFoundParID)
                curParToCompare = phoenixParameterList{curFoundParID(tmpPar,1),:};
                tmpParList{tmpPar,1} = strtrim(curParToCompare);
                if strcmp(tmpParList{tmpPar,1},curName)
                    IDs(curParameter,1) = curFoundParID(tmpPar,1);
                end
            end
            % Fallback
            if isnan(IDs(curParameter,1))
                IDs(curParameter,1) = curFoundParID(1);
                fprintf('Parameter %s found more than once and there is no exact match. Using first occurrence...\n', curName);
            end
        else
            % Parameter exists more than once
            IDs(curParameter,1) = NaN;
            fprintf('Parameter %s not found...\n', curName);
        end
        if ~isempty(IDs(curParameter,1)) && ~isnan(IDs(curParameter,1))
            parameters{newPar,2} = [phoenixParameterList{IDs(curParameter,1),2}];
            parameters{newPar,2} = strtrim(parameters{newPar,2});
            
            % Print out current information
            if debugMode
               fprintf('%s = %s\n',parameters{newPar,1},parameters{newPar,2}); 
            end
            % Iterate
            newPar = newPar+1;
        else
            % Parameter not found or multiple ones exists
            parameters{newPar,2} = '';
            % Iterate
            newPar = newPar+1;
        end
    end
end
