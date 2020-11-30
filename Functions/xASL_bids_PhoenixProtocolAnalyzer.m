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
    parameters = addParToList('tSequenceFileName',{},1);
    parameters = addParToList('UserScaleFactor',parameters,2);
    parameters = addParToList('M0Threshold',parameters,3);
    parameters = addParToList('TI1_us',parameters,4);
    parameters = addParToList('TI2_us',parameters,5);
    parameters = addParToList('sAsl.ulMode',parameters,6);
    parameters = addParToList('sAsl.ulAveragingMode',parameters,7);
    parameters = addParToList('sAsl.ulSuppressionMode',parameters,8);
    parameters = addParToList('sAsl.ulArrayLength',parameters,9);
    parameters = addParToList('M0Mode',parameters,10);
    parameters = addParToList('alTE[0]',parameters,11);
    parameters = addParToList('alTR[0]',parameters,12);
    parameters = addParToList('alTI[0]',parameters,13);
    parameters = addParToList('alTI[1]',parameters,14);
    parameters = addParToList('alTI[2]',parameters,15);
    parameters = addParToList('alTR.__attribute__.size',parameters,16);
    parameters = addParToList('alTI.__attribute__.size',parameters,17);
    parameters = addParToList('acFlowComp[0]',parameters,18);
    parameters = addParToList('lRepetitions',parameters,19);
    parameters = addParToList('lAverages',parameters,20);
    parameters = addParToList('lScanTimeSec',parameters,21);
    parameters = addParToList('CBFUpperLimit',parameters,22);
    parameters = addParToList('PerfPrepMode',parameters,23);
    parameters = addParToList('Slab_thickness_mm',parameters,24);
    parameters = addParToList('Dist_factor',parameters,25);
    parameters = addParToList('Flow_limit',parameters,26);
    parameters = addParToList('sGroupArray.sPSat.dGap',parameters,27);
    parameters = addParToList('sGroupArray.sPSat.dThickness',parameters,28);
    parameters = addParToList('sProtConsistencyInfo.tBaselineString',parameters,29);
    parameters = addParToList('sAsl.fFlowLimit',parameters,30);
    parameters = addParToList('sWipMemBlock.tFree',parameters,31);
    parameters = addParToList('sWipMemBlock.alFree.__attribute__.size',parameters,32);
    parameters = addParToList('sWipMemBlock.alFree[0]',parameters,33);
    parameters = addParToList('sWipMemBlock.alFree[1]',parameters,34);
    parameters = addParToList('sWipMemBlock.alFree[2]',parameters,35);
    parameters = addParToList('sWipMemBlock.alFree[5]',parameters,36);
    parameters = addParToList('sWipMemBlock.alFree[6]',parameters,37);
    parameters = addParToList('sWipMemBlock.alFree[63]',parameters,38);
    parameters = addParToList('sWipMemBlock.adFree.__attribute__.size',parameters,39);
    parameters = addParToList('sWipMemBlock.adFree[7]',parameters,40);
    parameters = addParToList('sWipMemBlock.adFree[8]',parameters,41);
    parameters = addParToList('sWipMemBlock.adFree[9]',parameters,42);
    parameters = addParToList('sWipMemBlock.adFree[10]',parameters,43);
    parameters = addParToList('sWipMemBlock.adFree[13]',parameters,44);

    
    %% Get the predefined parameters
    parameters = getPhoenixParameters(parameters,parameterList,false);
    
    %% Get xASL parameters
    parameters = convertCellArrayToStruct(parameters);

    %% Analyze parameters                                                   % Exemplary dataset in ExploreASL flavor library
    if contains(parameters.tSequenceFileName,'SiemensSeq')
        xasl.SequenceType = 'Siemens'; 
    elseif contains(parameters.tSequenceFileName,'CustomerSeq')
        xasl.SequenceType = 'Customer'; 
    end
    
    if contains(lower(parameters.tSequenceFileName),'ep2d')
        xasl.Sequence = '2D EPI'; 
    elseif contains(lower(parameters.tSequenceFileName),'gse') || ...
           contains(lower(parameters.tSequenceFileName),'grs3d')            % Siemens PASL 3DGRASE AD
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
    if contains(lower(parameters.sWipMemBlocktFree),'pcasl')                % Siemens PCASL 3DGRASE volunteer2
        xasl.LabelingType = 'PCASL';
    end
    
    if contains(lower(parameters.tSequenceFileName),'M0_') || ...           % Siemens PASL 3DGRASE AD
       contains(lower(parameters.tSequenceFileName),'_fid')                 % Siemens PASL 2D EPI noBsup AD
        xasl.ScanType = 'M0';
    end
    
    if parameters.sAslulMode==4
        xasl.Technique = 'Q2TIPS';
    end
    
    if ~isempty(parameters.alTI0) && ~isempty(parameters.alTI2)
        if parameters.alTI0~=100000 && parameters.alTI2~=7000000
            xasl.ScanType = 'ASL';
        elseif parameters.alTI0==100000 && parameters.alTI2==7000000
            xasl.ScanType = 'PseudoM0';
        end
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

% Add parameter to list
function parameters = addParToList(parName,parameters,parNum)
    parameters{parNum,1} = parName;
    parameters{parNum,2} = NaN;
end

