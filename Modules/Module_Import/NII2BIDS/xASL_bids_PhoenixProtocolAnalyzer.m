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
% 1. Extract parameters from Phoenix
% 2. Obtain basic information about the sequence
% 3. Reading sequence VEPCASL
% 4. Reading a default sequence
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          [bidsPar,sourcePar] = xASL_bids_PhoenixProtocolAnalyzer(parameterList);
%
% REFERENCES:       ...
% __________________________________
% Copyright @ 2015-2024 ExploreASL

    %% 1. Extract parameters from Phoenix
	parIndex = 1;
    sourcePar = addParToList('tSequenceFileName',{}, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('UserScaleFactor',sourcePar,parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('M0Threshold',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('TI1_us',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('TI2_us',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sAsl.ulMode',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sAsl.ulAveragingMode',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sAsl.ulSuppressionMode',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sAsl.ulArrayLength',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('M0Mode',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('alTE[0]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('alTR[0]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('alTI[0]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('alTI[1]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('alTI[2]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('alTR.__attribute__.size',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('alTI.__attribute__.size',sourcePar,parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('acFlowComp[0]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('lRepetitions',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('lAverages',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('lScanTimeSec',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('CBFUpperLimit',sourcePar,parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('PerfPrepMode',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('Slab_thickness_mm',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('Dist_factor',sourcePar,parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('Flow_limit',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sGroupArray.sPSat.dGap',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sGroupArray.sPSat.dThickness',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sProtConsistencyInfo.tBaselineString',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sAsl.fFlowLimit',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.tFree',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.alFree.__attribute__.size',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.alFree[0]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.alFree[1]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.alFree[2]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[3]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[4]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.alFree[5]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.alFree[6]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[9]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[10]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[11]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.alFree[63]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.adFree.__attribute__.size',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.adFree[0]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.adFree[1]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.adFree[7]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.adFree[8]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.adFree[9]',sourcePar,parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.adFree[10]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.adFree[13]',sourcePar, parIndex);parIndex = parIndex+1;

    % Get the predefined parameters
    sourcePar = getPhoenixParameters(sourcePar,parameterList,false);
    
    % Get xASL parameters
    sourcePar = convertCellArrayToStruct(sourcePar);

	%% 2. Obtain basic information about the sequence
	bidsPar = struct();
    % Find out if it's a product sequence or a patch
	if ~isempty(regexp(sourcePar.tSequenceFileName,'SiemensSeq', 'once'))
		sequenceType = 'Siemens'; % Not really a BIDS field
	elseif ~isempty(regexp(sourcePar.tSequenceFileName,'CustomerSeq', 'once'))
		sequenceType = 'Customer';
	end
    
	% Readout 2D or 3D
    if ~isempty(regexpi(sourcePar.tSequenceFileName, 'ep2d', 'once'))
        bidsPar.PulseSequenceType = '2D_EPI'; 
    elseif ~isempty(regexpi(sourcePar.tSequenceFileName, '(gse|grs3d', 'once'))
        bidsPar.PulseSequenceType = '3D_GRASE'; 
    end
    
	% Labeling type
    if ~isempty(regexpi(sourcePar.tSequenceFileName, 'pasl', 'once'))
        bidsPar.ArterialSpinLabelingType = 'PASL';
    end
    if ~isempty(regexpi(sourcePar.tSequenceFileName, 'casl', 'once'))
        bidsPar.ArterialSpinLabelingType = 'CASL';
    end
    if ~isempty(regexpi(sourcePar.tSequenceFileName, 'pcasl', 'once'))
        bidsPar.ArterialSpinLabelingType = 'PCASL';
    end
    if ~isempty(regexpi(sourcePar.sWipMemBlocktFree, 'pcasl', 'once'))                % Siemens PCASL 3DGRASE volunteer2
        bidsPar.ArterialSpinLabelingType = 'PCASL';
    end
    
	if ~isempty(regexpi(sourcePar.tSequenceFileName, '(M0_|_fid', 'once'))           % Siemens PASL 3DGRASE AD, Siemens PASL 2D EPI noBsup AD
        bidsPar.ScanType = 'M0';
	end
    
	bSequenceIdentified = false; % Try to identify different sequences and in case this doesn't work we rollback to the defaults

	%% 3. Reading sequence VEPCASL
	if ~bSequenceIdentified && ~isempty(regexp(sourcePar.tSequenceFileName,'VEPCASL', 'once'))
		if ~isempty(regexp(sourcePar.tSequenceFileName,'ep2d_VEPCASL', 'once'))
			% Basic check
			if ~isempty(sourcePar.sWipMemBlockalFree9) &&  ~isempty(sourcePar.sWipMemBlockalFree11) 
				% All basic passes checked, we can proceed to parameter parsing
				bSequenceIdentified = true;

				% Load the required parameters
				bidsPar.PostLabelingDelay = sourcePar.sWipMemBlockalFree11 / 1000;
				bidsPar.LabelingDuration = sourcePar.sWipMemBlockalFree11 / 1000;

				if ~isempty(sourcePar.sWipMemBlockalFree3) && sourcePar.sWipMemBlockalFree3 == 2
					bidsPar.BackgroundSuppression = true;
				end

				if ~isempty(sourcePar.sWipMemBlockalFree4)
					bidsPar.LabelingPulseFlipAngle = sourcePar.sWipMemBlockalFree4;
				end

				if ~isempty(sourcePar.sWipMemBlockalFree5)
					bidsPar.LabelingPulseDuration = sourcePar.sWipMemBlockalFree5 / 1000;
					if ~isempty(sourcePar.sWipMemBlockalFree6)
						bidsPar.LabelingPulseInterval = sourcePar.sWipMemBlockalFree6 / 1000 + bidsPar.LabelingPulseDuration;
					end
				end

				if ~isempty(sourcePar.sWipMemBlockadFree0)
					bidsPar.LabelingPulseAverageGradient = sourcePar.sWipMemBlockadFree0;
				end

				if ~isempty(sourcePar.sWipMemBlockadFree1)
					bidsPar.LabelingPulseMaximumGradient = sourcePar.sWipMemBlockadFree1;
				end
			else
				warning('Missing important parameters in Phoenix for Siemens VEPCASL sequence');
			end
		else
			% An unknown version detected, skipping to the default
			warning('An unknown version of Siemens VEPCASL sequence detected');
		end
	end

	%% 4. Reading a default sequence
	if ~bSequenceIdentified
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
				(~isempty(regexpi(bidsPar.SoftwareVersions,'N4_VB15A', 'once'))||...
				~isempty(regexpi(bidsPar.SoftwareVersions,'N4_VB17A', 'once'))||...
				~isempty(regexpi(bidsPar.SoftwareVersions,'N4_VB19A', 'once')))
			if isfield(sourcePar,'sGroupArraysPSatdGap') && ~isempty(sourcePar.sGroupArraysPSatdGap)
				bidsPar.LabelingDistance = sourcePar.sGroupArraysPSatdGap;
			end

			if isfield(bidsPar,'ArterialSpinLabelingType') && ~isempty(regexpi(bidsPar.ArterialSpinLabelingType,'pasl', 'once')) &&...
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
			if ~isempty(regexpi(bidsPar.ArterialSpinLabelingType,'pasl', 'once'))
				if isfield(bidsPar,'SoftwareVersions') &&...
						(~isempty(regexpi(bidsPar.SoftwareVersions,'N4_VB15A', 'once'))||...
						~isempty(regexpi(bidsPar.SoftwareVersions,'N4_VB17A', 'once'))||...
						~isempty(regexpi(bidsPar.SoftwareVersions,'N4_VB19A', 'once'))||...
						~isempty(regexpi(bidsPar.SoftwareVersions,'N4_VD13D', 'once'))||...
						~isempty(regexpi(bidsPar.SoftwareVersions,'N4_VE11C', 'once')))

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
			elseif ~isempty(regexpi(bidsPar.ArterialSpinLabelingType,'pcasl', 'once'))
				if isfield(sourcePar,'alTI0') && ~isempty(sourcePar.alTI0)
					if ~isfield(sourcePar, 'tSequenceFileName') || isempty(regexp(sourcePar.tSequenceFileName, 'ep2d_VEPCASL', 'once'))
						bidsPar.LabelingDuration = sourcePar.alTI0 / 1000 / 1000;
						if isfield(sourcePar,'alTI2') && ~isempty(sourcePar.alTI2)
							bidsPar.PostLabelingDelay = (sourcePar.alTI2-sourcePar.alTI0) / 1000 / 1000;
						end
					end
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

