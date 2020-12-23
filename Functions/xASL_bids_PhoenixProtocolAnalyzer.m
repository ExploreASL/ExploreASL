function [bidsPar,rawPar] = xASL_bids_PhoenixProtocolAnalyzer(parameterList)
%xASL_bids_PhoenixProtocolAnalyzer This function analyzes the parameter list of the phoenix protocol (tag = [0x29,0x1020]).
%
% FORMAT: [bidsPar,rawPar] = xASL_bids_PhoenixProtocolAnalyzer(parameterList);
%
% INPUT:
%        parameterList      - list of parameters from the reduced phoenix protocol (can be generated using xASL_bids_PhoenixProtocolReader) (REQUIRED)
%
% OUTPUT:
%        bidsPar               - list of BIDS parameters
%        rawPar         - list of raw unprocessed Phoenix parameters you want to know
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      This function analyzes the parameter list of the phoenix protocol (tag = [0x29,0x1020]).
%                   This function is usually called from xASL_bids_GetPhoenixProtocol.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          [bidsPar,rawPar] = xASL_bids_PhoenixProtocolAnalyzer(parameterList);
%
% REFERENCES:       ...
% __________________________________
% Copyright @ 2015-2020 ExploreASL

    %% Defaults
    rawPar = addParToList('tSequenceFileName',{},1);
    rawPar = addParToList('UserScaleFactor',rawPar,2);
    rawPar = addParToList('M0Threshold',rawPar,3);
    rawPar = addParToList('TI1_us',rawPar,4);
    rawPar = addParToList('TI2_us',rawPar,5);
    rawPar = addParToList('sAsl.ulMode',rawPar,6);
    rawPar = addParToList('sAsl.ulAveragingMode',rawPar,7);
    rawPar = addParToList('sAsl.ulSuppressionMode',rawPar,8);
    rawPar = addParToList('sAsl.ulArrayLength',rawPar,9);
    rawPar = addParToList('M0Mode',rawPar,10);
    rawPar = addParToList('alTE[0]',rawPar,11);
    rawPar = addParToList('alTR[0]',rawPar,12);
    rawPar = addParToList('alTI[0]',rawPar,13);
    rawPar = addParToList('alTI[1]',rawPar,14);
    rawPar = addParToList('alTI[2]',rawPar,15);
    rawPar = addParToList('alTR.__attribute__.size',rawPar,16);
    rawPar = addParToList('alTI.__attribute__.size',rawPar,17);
    rawPar = addParToList('acFlowComp[0]',rawPar,18);
    rawPar = addParToList('lRepetitions',rawPar,19);
    rawPar = addParToList('lAverages',rawPar,20);
    rawPar = addParToList('lScanTimeSec',rawPar,21);
    rawPar = addParToList('CBFUpperLimit',rawPar,22);
    rawPar = addParToList('PerfPrepMode',rawPar,23);
    rawPar = addParToList('Slab_thickness_mm',rawPar,24);
    rawPar = addParToList('Dist_factor',rawPar,25);
    rawPar = addParToList('Flow_limit',rawPar,26);
    rawPar = addParToList('sGroupArray.sPSat.dGap',rawPar,27);
    rawPar = addParToList('sGroupArray.sPSat.dThickness',rawPar,28);
    rawPar = addParToList('sProtConsistencyInfo.tBaselineString',rawPar,29);
    rawPar = addParToList('sAsl.fFlowLimit',rawPar,30);
    rawPar = addParToList('sWipMemBlock.tFree',rawPar,31);
    rawPar = addParToList('sWipMemBlock.alFree.__attribute__.size',rawPar,32);
    rawPar = addParToList('sWipMemBlock.alFree[0]',rawPar,33);
    rawPar = addParToList('sWipMemBlock.alFree[1]',rawPar,34);
    rawPar = addParToList('sWipMemBlock.alFree[2]',rawPar,35);
    rawPar = addParToList('sWipMemBlock.alFree[5]',rawPar,36);
    rawPar = addParToList('sWipMemBlock.alFree[6]',rawPar,37);
    rawPar = addParToList('sWipMemBlock.alFree[63]',rawPar,38);
    rawPar = addParToList('sWipMemBlock.adFree.__attribute__.size',rawPar,39);
    rawPar = addParToList('sWipMemBlock.adFree[7]',rawPar,40);
    rawPar = addParToList('sWipMemBlock.adFree[8]',rawPar,41);
    rawPar = addParToList('sWipMemBlock.adFree[9]',rawPar,42);
    rawPar = addParToList('sWipMemBlock.adFree[10]',rawPar,43);
    rawPar = addParToList('sWipMemBlock.adFree[13]',rawPar,44);

    
    %% Get the predefined parameters
    rawPar = getPhoenixParameters(rawPar,parameterList,false);
    
    %% Get xASL parameters
    rawPar = convertCellArrayToStruct(rawPar);

    %% Analyze parameters                                                   % Exemplary dataset in ExploreASL flavor library
    if contains(rawPar.tSequenceFileName,'SiemensSeq')
        bidsPar.SequenceType = 'Siemens'; % Not really a BIDS field
    elseif contains(rawPar.tSequenceFileName,'CustomerSeq')
        bidsPar.SequenceType = 'Customer'; 
	end
    
    if contains(lower(rawPar.tSequenceFileName),'ep2d')
        bidsPar.PulseSequenceType = '2D_EPI'; 
    elseif contains(lower(rawPar.tSequenceFileName),'gse') || ...
           contains(lower(rawPar.tSequenceFileName),'grs3d')            % Siemens PASL 3DGRASE AD
        bidsPar.PulseSequenceType = '3D_GRASE'; 
    end
    
    if contains(lower(rawPar.tSequenceFileName),'pasl')
        bidsPar.ArterialSpinLabelingType = 'PASL';
    end
    if contains(lower(rawPar.tSequenceFileName),'casl')
        bidsPar.ArterialSpinLabelingType = 'CASL';
    end
    if contains(lower(rawPar.tSequenceFileName),'pcasl')
        bidsPar.ArterialSpinLabelingType = 'PCASL';
    end
    if contains(lower(rawPar.sWipMemBlocktFree),'pcasl')                % Siemens PCASL 3DGRASE volunteer2
        bidsPar.ArterialSpinLabelingType = 'PCASL';
    end
    
    if contains(lower(rawPar.tSequenceFileName),'M0_') || ...           % Siemens PASL 3DGRASE AD
       contains(lower(rawPar.tSequenceFileName),'_fid')                 % Siemens PASL 2D EPI noBsup AD
        bidsPar.ScanType = 'M0';
	end
    
	%TODO
	% Saw this option also in QUIPSSII
    %if rawPar.sAslulMode==4
    %    bidsPar.BolusCutOffFlag = true;
	%	bidsPar.BolusCutOffTechnique = 'Q2TIPS';
    %end
    
	%TODO - does it really work for other sequences?
    if ~isempty(rawPar.alTI0) && ~isempty(rawPar.alTI2)
        if rawPar.alTI0~=100000 && rawPar.alTI2~=7000000
            bidsPar.ScanType = 'ASL';
        elseif rawPar.alTI0==100000 && rawPar.alTI2==7000000
            bidsPar.ScanType = 'PseudoM0';
        end
	end
    
	if isfield(rawPar,'sProtConsistencyInfotBaselineString') && ~isempty(rawPar.sProtConsistencyInfotBaselineString)
		bidsPar.SoftwareVersions = strrep(rawPar.sProtConsistencyInfotBaselineString,'"','');
	end
	
	if isfield(rawPar,'sGroupArraysPSatdGap') && ~isempty(rawPar.sGroupArraysPSatdGap)
		bidsPar.LabelingDistance = rawPar.sGroupArraysPSatdGap;
	end
	
	if isfield(bidsPar,'ArterialSpinLabelingType') && ~isempty(regexpi(bidsPar.ArterialSpinLabelingType,'pasl')) &&...
			isfield(rawPar,'sGroupArraysPSatdThickness') && ~isempty(rawPar.sGroupArraysPSatdThickness)
		bidsPar.LabelingSlabThickness = rawPar.sGroupArraysPSatdThickness;
	end
	
	% The function of this parameter is unclear
	%if isfield(rawPar,'Slab_thickness_mm') && ~isempty(rawPar.Slab_thickness_mm)
	%	bidsPar.LabelingSlabThickness = rawPar.Slab_thickness_mm;
	%end
	
	if isfield(rawPar,'alTE0') && ~isempty(rawPar.alTE0)
		bidsPar.EchoTime = rawPar.alTE0 / 1000 / 1000;
	end
	
	% If the labeling type is recognized, then proceed to labeling timing information extraction
	if isfield(bidsPar,'ArterialSpinLabelingType') 
		if~isempty(regexpi(bidsPar.ArterialSpinLabelingType,'pasl'))
			if isfield(bidsPar,'SoftwareVersions') &&...
					(~isempty(regexpi(bidsPar.SoftwareVersions,'N4_VB15A'))||...
					 ~isempty(regexpi(bidsPar.SoftwareVersions,'N4_VB17A'))||...
					 ~isempty(regexpi(bidsPar.SoftwareVersions,'N4_VB19A'))||...
					 ~isempty(regexpi(bidsPar.SoftwareVersions,'N4_VD13D'))||...
					 ~isempty(regexpi(bidsPar.SoftwareVersions,'N4_VE11C')))
				
				if isfield(rawPar,'alTI0') && ~isempty(rawPar.alTI0)
					bidsPar.PostLabelingDelay = rawPar.alTI0 / 1000 / 1000;
				end
				if isfield(rawPar,'alTI1') && ~isempty(rawPar.alTI1)
					if isfield(rawPar,'alTI2') && ~isempty(rawPar.alTI2)
						%N4_VB15A, N4_VB17A, N4_VB19A
						bidsPar.BolusCutOffDelayTime = [rawPar.alTI1, rawPar.alTI2] / 1000 / 1000;
						bidsPar.BolusCutOffFlag = true;
						bidsPar.BolusCutOffTechnique = 'Q2TIPS';
						bidsPar.PostLabelingDelay = rawPar.alTI0 / 1000 / 1000;
					else
						bidsPar.PostLabelingDelay = rawPar.alTI0 / 1000 / 1000;
						bidsPar.LabelingDuration = rawPar.alTI1 / 1000 / 1000;
					end
				else
					%N4_VD13D, N4_VE11C
					if isfield(rawPar,'alTI2') && ~isempty(rawPar.alTI2)
						bidsPar.PostLabelingDelay = rawPar.alTI2 / 1000 / 1000;
						bidsPar.BolusCutOffDelayTime = rawPar.alTI0 / 1000 / 1000;
						bidsPar.BolusCutOffFlag = true;
						bidsPar.BolusCutOffTechnique = 'QUIPSSII';
					end
				end
			end
		elseif ~isempty(regexpi(bidsPar.ArterialSpinLabelingType,'pasl'))
			bidsPar.LabelingDuration = rawPar.alTI0 / 1000 / 1000;
		elseif ~isempty(regexpi(bidsPar.ArterialSpinLabelingType,'pcasl'))
			if isfield(rawPar,'alTI0') && ~isempty(rawPar.alTI0)
				bidsPar.LabelingDuration = rawPar.alTI0 / 1000 / 1000;
				if isfield(rawPar,'alTI2') && ~isempty(rawPar.alTI2)
					bidsPar.PostLabelingDelay = (rawPar.alTI2-rawPar.alTI0) / 1000 / 1000;
				end
			end
		end
	end

end

%% Get individual phoenix parameter
function value = getPhoePar(rawPar,curParToExtract)
    % Get cell ID
    cellID = find(contains(rawPar, curParToExtract));
    % Return value
    value = rawPar{cellID,2};
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
			if debugMode
				fprintf('Parameter %s not found...\n', curName);
			end
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

