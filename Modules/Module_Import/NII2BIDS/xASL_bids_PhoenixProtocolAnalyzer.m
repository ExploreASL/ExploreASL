function [bidsPar,sourcePar] = xASL_bids_PhoenixProtocolAnalyzer(parameterList)
%xASL_bids_PhoenixProtocolAnalyzer This function analyzes the parameter list of the phoenix protocol (tag = [0x29,0x1020]).
%
% FORMAT: [bidsPar,sourcePar] = xASL_bids_PhoenixProtocolAnalyzer(parameterList);
%
% INPUT:
%        parameterList      - list of parameters from the reduced phoenix protocol (can be generated using xASL_bids_PhoenixProtocolReader) (REQUIRED)
%
% OUTPUT:
%        bidsPar           - list of BIDS parameters
%        sourcePar         - list of raw unprocessed Phoenix parameters you want to know
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      This function analyzes the parameter list of the phoenix protocol (tag = [0x29,0x1020]).
%                   This function is usually called from xASL_bids_GetPhoenixProtocol.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          [bidsPar,sourcePar] = xASL_bids_PhoenixProtocolAnalyzer(parameterList);
%
% REFERENCES:       ...
% __________________________________
% Copyright @ 2015-2021 ExploreASL

    %% Defaults
    sourcePar = addParToList('tSequenceFileName',{},1);
    sourcePar = addParToList('UserScaleFactor',sourcePar,2);
    sourcePar = addParToList('M0Threshold',sourcePar,3);
    sourcePar = addParToList('TI1_us',sourcePar,4);
    sourcePar = addParToList('TI2_us',sourcePar,5);
    sourcePar = addParToList('sAsl.ulMode',sourcePar,6);
    sourcePar = addParToList('sAsl.ulAveragingMode',sourcePar,7);
    sourcePar = addParToList('sAsl.ulSuppressionMode',sourcePar,8);
    sourcePar = addParToList('sAsl.ulArrayLength',sourcePar,9);
    sourcePar = addParToList('M0Mode',sourcePar,10);
    sourcePar = addParToList('alTE[0]',sourcePar,11);
    sourcePar = addParToList('alTR[0]',sourcePar,12);
    sourcePar = addParToList('alTI[0]',sourcePar,13);
    sourcePar = addParToList('alTI[1]',sourcePar,14);
    sourcePar = addParToList('alTI[2]',sourcePar,15);
    sourcePar = addParToList('alTR.__attribute__.size',sourcePar,16);
    sourcePar = addParToList('alTI.__attribute__.size',sourcePar,17);
    sourcePar = addParToList('acFlowComp[0]',sourcePar,18);
    sourcePar = addParToList('lRepetitions',sourcePar,19);
    sourcePar = addParToList('lAverages',sourcePar,20);
    sourcePar = addParToList('lScanTimeSec',sourcePar,21);
    sourcePar = addParToList('CBFUpperLimit',sourcePar,22);
    sourcePar = addParToList('PerfPrepMode',sourcePar,23);
    sourcePar = addParToList('Slab_thickness_mm',sourcePar,24);
    sourcePar = addParToList('Dist_factor',sourcePar,25);
    sourcePar = addParToList('Flow_limit',sourcePar,26);
    sourcePar = addParToList('sGroupArray.sPSat.dGap',sourcePar,27);
    sourcePar = addParToList('sGroupArray.sPSat.dThickness',sourcePar,28);
    sourcePar = addParToList('sProtConsistencyInfo.tBaselineString',sourcePar,29);
    sourcePar = addParToList('sAsl.fFlowLimit',sourcePar,30);
    sourcePar = addParToList('sWipMemBlock.tFree',sourcePar,31);
    sourcePar = addParToList('sWipMemBlock.alFree.__attribute__.size',sourcePar,32);
    sourcePar = addParToList('sWipMemBlock.alFree[0]',sourcePar,33);
    sourcePar = addParToList('sWipMemBlock.alFree[1]',sourcePar,34);
    sourcePar = addParToList('sWipMemBlock.alFree[2]',sourcePar,35);
    sourcePar = addParToList('sWipMemBlock.alFree[5]',sourcePar,36);
    sourcePar = addParToList('sWipMemBlock.alFree[6]',sourcePar,37);
    sourcePar = addParToList('sWipMemBlock.alFree[63]',sourcePar,38);
    sourcePar = addParToList('sWipMemBlock.adFree.__attribute__.size',sourcePar,39);
    sourcePar = addParToList('sWipMemBlock.adFree[7]',sourcePar,40);
    sourcePar = addParToList('sWipMemBlock.adFree[8]',sourcePar,41);
    sourcePar = addParToList('sWipMemBlock.adFree[9]',sourcePar,42);
    sourcePar = addParToList('sWipMemBlock.adFree[10]',sourcePar,43);
    sourcePar = addParToList('sWipMemBlock.adFree[13]',sourcePar,44);

    
    %% Get the predefined parameters
    sourcePar = getPhoenixParameters(sourcePar,parameterList,false);
    
    %% Get xASL parameters
    sourcePar = convertCellArrayToStruct(sourcePar);

	bidsPar = struct();
    %% Analyze parameters                                                   % Exemplary dataset in ExploreASL flavor library
	if ~isempty(strfind(sourcePar.tSequenceFileName,'SiemensSeq'))
		bidsPar.SequenceType = 'Siemens'; % Not really a BIDS field
	elseif ~isempty(strfind(sourcePar.tSequenceFileName,'CustomerSeq'))
		bidsPar.SequenceType = 'Customer';
	end
    
    if ~isempty(strfind(lower(sourcePar.tSequenceFileName),'ep2d'))
        bidsPar.PulseSequenceType = '2D_EPI'; 
    elseif ~isempty(strfind(lower(sourcePar.tSequenceFileName),'gse')) || ...
           ~isempty(strfind(lower(sourcePar.tSequenceFileName),'grs3d'))            % Siemens PASL 3DGRASE AD
        bidsPar.PulseSequenceType = '3D_GRASE'; 
    end
    
    if ~isempty(strfind(lower(sourcePar.tSequenceFileName),'pasl'))
        bidsPar.ArterialSpinLabelingType = 'PASL';
    end
    if ~isempty(strfind(lower(sourcePar.tSequenceFileName),'casl'))
        bidsPar.ArterialSpinLabelingType = 'CASL';
    end
    if ~isempty(strfind(lower(sourcePar.tSequenceFileName),'pcasl'))
        bidsPar.ArterialSpinLabelingType = 'PCASL';
    end
    if ~isempty(strfind(lower(sourcePar.sWipMemBlocktFree),'pcasl'))                % Siemens PCASL 3DGRASE volunteer2
        bidsPar.ArterialSpinLabelingType = 'PCASL';
    end
    
    if ~isempty(strfind(lower(sourcePar.tSequenceFileName),'M0_')) || ...           % Siemens PASL 3DGRASE AD
       ~isempty(strfind(lower(sourcePar.tSequenceFileName),'_fid'))                 % Siemens PASL 2D EPI noBsup AD
        bidsPar.ScanType = 'M0';
	end
    
	%TODO
	% Saw this option also in QUIPSSII
    %if sourcePar.sAslulMode==4
    %    bidsPar.BolusCutOffFlag = true;
	%	bidsPar.BolusCutOffTechnique = 'Q2TIPS';
    %end
    
	%TODO - does it really work for other sequences?
    if ~isempty(sourcePar.alTI0) && ~isempty(sourcePar.alTI2)
        if sourcePar.alTI0~=100000 && sourcePar.alTI2~=7000000
            bidsPar.ScanType = 'ASL';
        elseif sourcePar.alTI0==100000 && sourcePar.alTI2==7000000
            bidsPar.ScanType = 'PseudoM0';
        end
	end
    
	if isfield(sourcePar,'sProtConsistencyInfotBaselineString') && ~isempty(sourcePar.sProtConsistencyInfotBaselineString)
		bidsPar.SoftwareVersions = strrep(sourcePar.sProtConsistencyInfotBaselineString,'"','');
	end
	
	% N4_VD13D and N4_VE11C have also slab_thickness_mm parameter defined and it is unclear if this is imaging FOV or labeling slab thickness...
	if isfield(bidsPar,'SoftwareVersions') &&...
			(~isempty(regexpi(bidsPar.SoftwareVersions,'N4_VB15A'))||...
			 ~isempty(regexpi(bidsPar.SoftwareVersions,'N4_VB17A'))||...
			 ~isempty(regexpi(bidsPar.SoftwareVersions,'N4_VB19A')))
		if isfield(sourcePar,'sGroupArraysPSatdGap') && ~isempty(sourcePar.sGroupArraysPSatdGap)
			bidsPar.LabelingDistance = sourcePar.sGroupArraysPSatdGap;
		end
	
		if isfield(bidsPar,'ArterialSpinLabelingType') && ~isempty(regexpi(bidsPar.ArterialSpinLabelingType,'pasl')) &&...
				isfield(sourcePar,'sGroupArraysPSatdThickness') && ~isempty(sourcePar.sGroupArraysPSatdThickness)
			bidsPar.LabelingSlabThickness = sourcePar.sGroupArraysPSatdThickness;
		end
	end
	% The function of this parameter is unclear
	%if isfield(sourcePar,'Slab_thickness_mm') && ~isempty(sourcePar.Slab_thickness_mm)
	%	bidsPar.LabelingSlabThickness = sourcePar.Slab_thickness_mm;
	%end
	
	if isfield(sourcePar,'alTE0') && ~isempty(sourcePar.alTE0)
		bidsPar.EchoTime = sourcePar.alTE0 / 1000 / 1000;
		lastEchoTime = bidsPar.EchoTime;
		lastEchoTimeCount = 1;
		while lastEchoTime > 0
			lastEchoTimeCount = lastEchoTimeCount + 1;
			echoTimePar = getPhoenixParameters({['alTE[' num2str(lastEchoTimeCount) ']']},parameterList,0);
			if isempty(echoTimePar{2})
				lastEchoTime = 0;
			else
				lastEchoTime = str2num(echoTimePar{2})/1000/1000;
				bidsPar.EchoTime(lastEchoTimeCount) = lastEchoTime;
			end
		end		
		
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
				
				if isfield(sourcePar,'alTI0') && ~isempty(sourcePar.alTI0)
					bidsPar.PostLabelingDelay = sourcePar.alTI0 / 1000 / 1000;
				end
				if isfield(sourcePar,'alTI1') && ~isempty(sourcePar.alTI1)
					if isfield(sourcePar,'alTI2') && ~isempty(sourcePar.alTI2)
						%N4_VB15A, N4_VB17A, N4_VB19A
						bidsPar.BolusCutOffDelayTime = [sourcePar.alTI1, sourcePar.alTI2] / 1000 / 1000;
						bidsPar.BolusCutOffFlag = true;
						bidsPar.BolusCutOffTechnique = 'Q2TIPS';
						bidsPar.PostLabelingDelay = sourcePar.alTI0 / 1000 / 1000;
					else
						bidsPar.PostLabelingDelay = sourcePar.alTI0 / 1000 / 1000;
						bidsPar.LabelingDuration = sourcePar.alTI1 / 1000 / 1000;
					end
				else
					%N4_VD13D, N4_VE11C
					if isfield(sourcePar,'alTI2') && ~isempty(sourcePar.alTI2)
						bidsPar.PostLabelingDelay = sourcePar.alTI2 / 1000 / 1000;
						bidsPar.BolusCutOffDelayTime = sourcePar.alTI0 / 1000 / 1000;
						bidsPar.BolusCutOffFlag = true;
						bidsPar.BolusCutOffTechnique = 'QUIPSSII';
					end
				end
			end
		elseif ~isempty(regexpi(bidsPar.ArterialSpinLabelingType,'pasl'))
			bidsPar.LabelingDuration = sourcePar.alTI0 / 1000 / 1000;
		elseif ~isempty(regexpi(bidsPar.ArterialSpinLabelingType,'pcasl'))
			if isfield(sourcePar,'alTI0') && ~isempty(sourcePar.alTI0)
				bidsPar.LabelingDuration = sourcePar.alTI0 / 1000 / 1000;
				if isfield(sourcePar,'alTI2') && ~isempty(sourcePar.alTI2)
					bidsPar.PostLabelingDelay = (sourcePar.alTI2-sourcePar.alTI0) / 1000 / 1000;
				end
			end
		end
	end

end

%% Get individual phoenix parameter
function value = getPhoePar(sourcePar,curParToExtract)
    % Get cell ID
    IndexC = strfind(sourcePar,curParToExtract);
    cellID = find(not(cellfun('isempty',IndexC)));
    % Return value
    value = sourcePar{cellID,2};
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
        IndexC = strfind(phoenixParameterList,curName);
        curFoundParID = find(not(cellfun('isempty',IndexC)));
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

