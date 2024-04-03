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
% 4. Reading the _VE11C sequences from DJJ Wang
% 5. Reading the CSL_818 sequences
% 6. Reading a default sequence
% 7. Reading TE vector and TR vector
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          [bidsPar,sourcePar] = xASL_bids_PhoenixProtocolAnalyzer(parameterList);
%
% __________________________________
% Copyright @ 2015-2024 ExploreASL

    %% 1. Extract parameters from Phoenix - create a list of parameters
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
	sourcePar = addParToList('sFastImaging.lSegments',sourcePar, parIndex);parIndex = parIndex+1;
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
	sourcePar = addParToList('sWipMemBlock.alFree[7]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[8]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[9]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[10]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[11]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[12]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[13]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[14]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[15]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[16]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[17]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[18]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[19]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[20]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[30]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[31]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[32]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[33]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[34]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[35]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[36]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[37]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[38]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[39]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.alFree[40]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.alFree[63]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.adFree.__attribute__.size',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.adFree[0]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.adFree[1]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.adFree[2]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.adFree[3]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.adFree[4]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.adFree[5]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.adFree[6]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.adFree[7]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.adFree[8]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.adFree[9]',sourcePar,parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.adFree[10]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.adFree[11]',sourcePar, parIndex);parIndex = parIndex+1;
	sourcePar = addParToList('sWipMemBlock.adFree[12]',sourcePar, parIndex);parIndex = parIndex+1;
    sourcePar = addParToList('sWipMemBlock.adFree[13]',sourcePar, parIndex);parIndex = parIndex+1;

    % Get the predefined parameters
    sourcePar = xASL_bids_PhoenixProtocolAnalyzer_getPhoenixParameters(sourcePar,parameterList,false);
    
    % Get xASL parameters
    sourcePar = xASL_bids_PhoenixProtocolAnalyzer_convertCellArrayToStruct(sourcePar);

	%% 2. Obtain basic information about the sequence
	bidsPar = struct();
    % Find out if it's a product sequence or a patch
	if ~isempty(regexpi(sourcePar.tSequenceFileName,'SiemensSeq', 'once'))
		sequenceType = 'Siemens'; % Not really a BIDS field
	elseif ~isempty(regexpi(sourcePar.tSequenceFileName,'CustomerSeq', 'once'))
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
    
	% Read the software version
	if isfield(sourcePar,'sProtConsistencyInfotBaselineString') && ~isempty(sourcePar.sProtConsistencyInfotBaselineString)
		bidsPar.SoftwareVersions = strrep(sourcePar.sProtConsistencyInfotBaselineString,'"','');
	end

	bSequenceIdentified = false; % Try to identify different sequences and in case this doesn't work we rollback to the defaults

	%% 3. Reading sequence VEPCASL
	if ~bSequenceIdentified && ~isempty(regexp(sourcePar.tSequenceFileName,'VEPCASL', 'once'))
		% Check that the PCASL tag is filled as ON (==1)
		if isempty(sourcePar.sWipMemBlockalFree0) || sourcePar.sWipMemBlockalFree0 ~= 1
			error('VEPCASL sequence without PCASL activated');
		end
		% 2D-EPI VEPCASL
		if ~isempty(regexp(sourcePar.tSequenceFileName,'to_ep2d_VEPCASL', 'once'))
			% Check that PLD and LabelingDuration tags are provided
			if ~isempty(sourcePar.sWipMemBlockalFree9) &&  ~isempty(sourcePar.sWipMemBlockalFree11) 
				% Proceed to parameter parsing

				% Load the required parameters
				bidsPar.LabelingDuration = sourcePar.sWipMemBlockalFree9 / 1000;

				% PLDs are in alFree11 till alFree20
				bidsPar.PostLabelingDelay = sourcePar.sWipMemBlockalFree11 / 1000;
				iPLD = 2;
				while iPLD <= 10 && ~isempty(sourcePar.(['sWipMemBlockalFree' num2str(iPLD+10)]))
					bidsPar.PostLabelingDelay(iPLD) = sourcePar.(['sWipMemBlockalFree' num2str(iPLD+10)]) / 1000;
					iPLD = iPLD + 1;
				end

				% Check if all control and label images are acquired and saved
				if sourcePar.sWipMemBlockalFree2 == 4
					% If more than 1 PLD is given, we have to double that for control and labels
					if length(bidsPar.PostLabelingDelay) > 1
						bidsPar.PostLabelingDelay = [bidsPar.PostLabelingDelay;bidsPar.PostLabelingDelay];
						bidsPar.PostLabelingDelay = bidsPar.PostLabelingDelay(:)';
					end
				end
				
				% Background suppression
				if ~isempty(sourcePar.sWipMemBlockalFree3) && (sourcePar.sWipMemBlockalFree3 == 2 || sourcePar.sWipMemBlockalFree3 == 3)
					bidsPar.BackgroundSuppression = true;
				end

				% Flip angle
				if ~isempty(sourcePar.sWipMemBlockalFree4)
					bidsPar.LabelingPulseFlipAngle = sourcePar.sWipMemBlockalFree4;
				end

				% Details of the labeling pulse design
				if ~isempty(sourcePar.sWipMemBlockalFree5)
					bidsPar.LabelingPulseDuration = sourcePar.sWipMemBlockalFree5 / 1000;
					if ~isempty(sourcePar.sWipMemBlockalFree6)
						bidsPar.LabelingPulseInterval = sourcePar.sWipMemBlockalFree6 / 1000 + bidsPar.LabelingPulseDuration;
					end
				end
			else
				warning('Missing important parameters in Phoenix for Siemens VEPCASL sequence');
			end

		elseif ~isempty(regexp(sourcePar.tSequenceFileName, 'jw_tgse_VEPCASL', 'once'))
			% Check that PLD and LabelingDuration tags are provided
			if ~isempty(sourcePar.sWipMemBlockalFree10) &&  ~isempty(sourcePar.sWipMemBlockalFree31) 

				% Load the required parameters
				bidsPar.LabelingDuration = sourcePar.sWipMemBlockalFree10 / 1000;

				% For 3D GRASE, PLDs are in alFree31 till alFree40
				bidsPar.PostLabelingDelay = sourcePar.sWipMemBlockalFree30 / 1000;
				iPLD = 2;
				while iPLD <= 10 && ~isempty(sourcePar.(['sWipMemBlockalFree' num2str(iPLD+29)]))
					bidsPar.PostLabelingDelay(iPLD) = sourcePar.(['sWipMemBlockalFree' num2str(iPLD+29)]) / 1000;
					iPLD = iPLD + 1;
				end

				% Check if all control and label images are acquired and saved
				if sourcePar.sWipMemBlockalFree4 == 4
					% If more than 1 PLD is given, we have to double that for control and labels
					if length(bidsPar.PostLabelingDelay) > 1
						bidsPar.PostLabelingDelay = [bidsPar.PostLabelingDelay;bidsPar.PostLabelingDelay];
						bidsPar.PostLabelingDelay = bidsPar.PostLabelingDelay(:)';
					end
				end
				
				% Background suppression
				if ~isempty(sourcePar.sWipMemBlockalFree5) && (sourcePar.sWipMemBlockalFree5 == 3 || sourcePar.sWipMemBlockalFree5 == 2)
					bidsPar.BackgroundSuppression = true;
				end

				% Flip angle
				if ~isempty(sourcePar.sWipMemBlockalFree6)
					bidsPar.LabelingPulseFlipAngle = sourcePar.sWipMemBlockalFree6;
				end

				% Details of the labeling pulse design
				if ~isempty(sourcePar.sWipMemBlockalFree7)
					bidsPar.LabelingPulseDuration = sourcePar.sWipMemBlockalFree7/ 1000;
					if ~isempty(sourcePar.sWipMemBlockalFree8)
						bidsPar.LabelingPulseInterval = sourcePar.sWipMemBlockalFree8 / 1000 + bidsPar.LabelingPulseDuration;
					end
				end
			else
				warning('Missing important parameters in Phoenix for Siemens VEPCASL sequence');
			end	
		else
			% An unknown version detected, skipping to the default
			error('An unknown version of Siemens VEPCASL sequence detected');
		end

		bSequenceIdentified = true;

		% Common parameters for all VEPCASL implementations
		if ~isempty(sourcePar.sWipMemBlockalFree0) && sourcePar.sWipMemBlockalFree0 == 1
			bidsPar.ArterialSpinLabelingType = 'PCASL';
		end

		if ~isempty(sourcePar.sWipMemBlockadFree0)
			bidsPar.LabelingPulseAverageGradient = sourcePar.sWipMemBlockadFree0;
		end

		if ~isempty(sourcePar.sWipMemBlockadFree1)
			bidsPar.LabelingPulseMaximumGradient = sourcePar.sWipMemBlockadFree1;
		end
	end

	%% 4. Reading the _VE11C sequences from DJJ Wang
	if ~bSequenceIdentified && ~isempty(regexpi(sourcePar.tSequenceFileName,'pcasl_ve11c', 'once'))

		if ~isempty(regexpi(sourcePar.tSequenceFileName,'ep2d_pcasl_ve11c', 'once'))
			% 2DEPI VE11C PCASL
			if ~isempty(sourcePar.sWipMemBlockadFree2)
				bidsPar.PostLabelingDelay = sourcePar.sWipMemBlockadFree2 / 1000000.0;
			end

		elseif ~isempty(regexpi(sourcePar.tSequenceFileName,'tgse_pcasl_ve11c', 'once'))
			% 3DGRASE VE11C PCASL
			if ~isempty(sourcePar.sWipMemBlockadFree2)
				bidsPar.PostLabelingDelay = sourcePar.sWipMemBlockadFree2 / 1000000.0;
			end

			if ~isempty(sourcePar.sWipMemBlockalFree13) && sourcePar.sWipMemBlockalFree13 == 1
				bidsPar.BackgroundSuppression = true;
			else
				bidsPar.BackgroundSuppression = false;
			end

		else
			error('Unknown variant of PCASL_VE11C');
		end

		if ~isempty(sourcePar.sWipMemBlockadFree1)
			bidsPar.LabelingDistance = sourcePar.sWipMemBlockadFree1;
		end

		bSequenceIdentified = true;
		bidsPar.ArterialSpinLabelingType = 'PCASL';

		if ~isempty(sourcePar.sWipMemBlockadFree10)
			bidsPar.LabelingPulseAverageGradient = sourcePar.sWipMemBlockadFree10 / 10.0; % MeanGz x 10 mT/m
		end

		%sourcePar.sWipMemBlockadFree4 is RFGap in micro seconds (typically 360)
		%sourcePar.sWipMemBlockadFree11 is PhiAdjust 100
		%sourcePar.sWipMemBlockadFree11 T1 in usec

		if ~isempty(sourcePar.sWipMemBlockadFree3)
			%  Number of 20 RF pulses (each block is 18.4ms, recommend 82 for total labeling duration of 1.5s
			bidsPar.LabelingDuration = sourcePar.sWipMemBlockadFree3 * 18.4 / 1000.0;
		end
		
		% M0 TR for this sequence is 2.0s by default
		bidsPar.RepetitionTimePreparationM0 = 2.0;
	end
	
	%% 5. Reading the CSL_818 sequences
	% Identify the sequence
	if ~bSequenceIdentified && ~isempty(regexpi(sourcePar.tSequenceFileName, 'csl_818', 'once'))
		if ~bSequenceIdentified && ~isempty(regexpi(sourcePar.tSequenceFileName, 'csl_818', 'once'))
			% Read Flip angle of the labeling pulse
			if ~isempty(sourcePar.sWipMemBlockadFree13)
				bidsPar.LabelingPulseFlipAngle = sourcePar.sWipMemBlockadFree13;
			end

			% Read PLD and LabDur from the common Phoenix fields
			if ~isempty(sourcePar.alTI0)
				bidsPar.LabelingDuration = sourcePar.alTI0 / 1000 / 1000;
				if ~isempty(sourcePar.alTI2)
					bidsPar.PostLabelingDelay = (sourcePar.alTI2-sourcePar.alTI0) / 1000 / 1000;
				end
			end

			% Number of control/label pairs
			if ~isempty(sourcePar.lRepetitions)
				bidsPar.NumberOfAverages = sourcePar.lRepetitions + 1;
			end

			% Number of segments in the 3D GRASE readout
			if ~isempty(sourcePar.sFastImaginglSegments)
				bidsPar.NumberSegments = sourcePar.sFastImaginglSegments;
			end
 
			% The repetitiontime of M0
			bidsPar.RepetitionTimePreparationM0 = sourcePar.sWipMemBlockalFree2 / 1000;

			% The normal repetition time
			bidsPar.RepetitionTimePreparation = sourcePar.alTR0 / 1000 / 1000;
		else
			error('Unknown variant of CSL_818');
		end
		bSequenceIdentified = true;
	end
	%% 6. Reading a default sequence
	if ~bSequenceIdentified
		% TODO
		% Saw this option also in QUIPSSII
		%if sourcePar.sAslulMode==4
		%    bidsPar.BolusCutOffFlag = true;
		%	bidsPar.BolusCutOffTechnique = 'Q2TIPS';
		%end

		% TODO - does it really work for other sequences?
		if ~isempty(sourcePar.alTI0) && ~isempty(sourcePar.alTI2)
			if sourcePar.alTI0~=100000 && sourcePar.alTI2~=7000000
				bidsPar.ScanType = 'ASL';
			elseif sourcePar.alTI0==100000 && sourcePar.alTI2==7000000
				bidsPar.ScanType = 'PseudoM0';
			end
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
					bidsPar.LabelingDuration = sourcePar.alTI0 / 1000 / 1000;
					if isfield(sourcePar,'alTI2') && ~isempty(sourcePar.alTI2)
						bidsPar.PostLabelingDelay = (sourcePar.alTI2-sourcePar.alTI0) / 1000 / 1000;
					end
				end
			end
		end
	end

	%% 7. Reading TE vector and TR vector
	if isfield(sourcePar,'alTE0') && ~isempty(sourcePar.alTE0)
		bidsPar.EchoTime = sourcePar.alTE0 / 1000 / 1000;
		lastEchoTime = bidsPar.EchoTime;
		lastEchoTimeCount = 1;
		while lastEchoTime > 0
			lastEchoTimeCount = lastEchoTimeCount + 1;
			echoTimePar = xASL_bids_PhoenixProtocolAnalyzer_getPhoenixParameters({['alTE[' num2str(lastEchoTimeCount-1) ']']},parameterList,0);
			% If the fields are empty or start decreasing, then we stop reading
			if isempty(echoTimePar{2}) || (lastEchoTime > str2num(echoTimePar{2})/1000/1000)
				lastEchoTime = 0;
			else
				lastEchoTime = str2num(echoTimePar{2})/1000/1000;
				bidsPar.EchoTime(lastEchoTimeCount) = lastEchoTime;
			end
		end
	end

end
%% --------------------------------------------------------------------------------------------------------
%% Get individual phoenix parameter
%% --------------------------------------------------------------------------------------------------------
function value = xASL_bids_PhoenixProtocolAnalyzer_getPhoePar(sourcePar,curParToExtract)
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

%% --------------------------------------------------------------------------------------------------------
%% Convert cell array to struct
%% --------------------------------------------------------------------------------------------------------
function parameterStruct = xASL_bids_PhoenixProtocolAnalyzer_convertCellArrayToStruct(parameters)
    % Create struct
    parameterStruct = struct;
    % Get number of parameters
    parameterNames = parameters(:,1);
    numberOfParameters = numel(parameterNames);
    % Convert each line to a field
    for curParameter=1:numberOfParameters
        curParToExtract = parameters{curParameter,1};
        value = xASL_bids_PhoenixProtocolAnalyzer_getPhoePar(parameters,curParToExtract);
        validFieldName = matlab.lang.makeValidName(curParToExtract,'ReplacementStyle','delete');
        parameterStruct.(validFieldName) = value;        
    end
end

%% --------------------------------------------------------------------------------------------------------
%% Get ID of parameter name
%% --------------------------------------------------------------------------------------------------------
function parameters = xASL_bids_PhoenixProtocolAnalyzer_getPhoenixParameters(parameters,phoenixParameterList,debugMode)

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

