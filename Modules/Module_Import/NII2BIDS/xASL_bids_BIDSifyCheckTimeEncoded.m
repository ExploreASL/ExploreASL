function [jsonOut,bTimeEncoded, bTimeEncodedFME] = xASL_bids_BIDSifyCheckTimeEncoded(jsonIn, jsonOut, nVolumes)
%xASL_bids_BIDSifyCheckTimeEncoded Check for time encoded sequence
%
% FORMAT: [jsonOut,bTimeEncoded, bTimeEncodedFME] = xASL_bids_BIDSifyCheckTimeEncoded(jsonIn, jsonOut, nVolumes)
%
% INPUT:
%   jsonIn    - JSON with the input fields - from DICOMs (REQUIRED)
%   jsonOut   - Output JSON in progress from the parent function xASL_bids_BIDSifyASLJSON (REQUIRED)
%   nVolumes  - number of acquired volumes (REQUIRED)
%
% OUTPUT: 
%   jsonOut   - Output JSON in progress from the parent function xASL_bids_BIDSifyASLJSON
%   bTimeEncoded    - Boolean describing if the current sequence is a time encoded sequence
%   bTimeEncodedFME - Boolean describing if the current sequence is a specific FME time encoded sequence
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Check for time encoded sequence.
%
% EXAMPLE:     [jsonOut,bTimeEncoded,bTimeEncodedFME] = xASL_bids_BIDSifyCheckTimeEncoded(jsonIn,jsonOut);
%
% __________________________________
% Copyright 2015-2021 ExploreASL

    if nargin < 3
		error('Three input parameters required');
	end
	
    % Check for time encoded sequences
    if isfield(jsonOut,'TimeEncodedMatrixSize') && ~isempty(jsonOut.TimeEncodedMatrixSize) || ... % Should be 4, 8 or 12
        isfield(jsonOut,'TimeEncodedMatrixType') % Natural or walsh
        bTimeEncoded = true; 
    else 
        bTimeEncoded = false; 
    end

    % Check for specific time encoded sequence of FME (Fraunhofer Mevis)
	bTimeEncodedFME = xASL_imp_CheckIfFME(jsonIn, jsonOut, bTimeEncoded);
	
    if bTimeEncodedFME
        bTimeEncoded = true;    
    end
    
	if isfield(jsonOut, 'EchoTime')
		NumberEchoTimes = length(uniquetol(jsonOut.EchoTime,0.001)); % Obtain the number of echo times

		% Either 1 TE or matching the number of the volumes
		if length(jsonOut.EchoTime) > 1 && nVolumes>length(jsonOut.EchoTime)
			if mod(nVolumes,length(jsonOut.EchoTime)) == 0
				switch(jsonOut.Manufacturer)
					case {'Siemens'}
						jsonOut.EchoTime = repmat(jsonOut.EchoTime(:), 1, nVolumes/length(jsonOut.EchoTime))';
						jsonOut.EchoTime = jsonOut.EchoTime(:);
					case {'Philips'}
						jsonOut.EchoTime = repmat(jsonOut.EchoTime(:), nVolumes/length(jsonOut.EchoTime),1)';
						jsonOut.EchoTime = jsonOut.EchoTime(:);
					otherwise
						warning(['Cannot resort multi-TE for ' jsouOut.Manufacturer]);
				end
			else
				warning('Number of EchoTimes %d does not match the number of volumes %d\n',length(jsonOut.EchoTime),nVolumes);
			end
		end
	else
		NumberEchoTimes = 1;
	end
		
    if bTimeEncoded
        if isfield(jsonOut, 'PostLabelingDelay')
            % From the import, the length of EchoTime should correspond to the number of volumes
            if length(jsonOut.PostLabelingDelay) > 1 && length(jsonOut.PostLabelingDelay) < nVolumes
                % So here, we first make sure that each PLD is repeated for the whole block of echo-times
				if mod(nVolumes,length(jsonOut.PostLabelingDelay)*NumberEchoTimes) == 0
					switch (jsonOut.Manufacturer)
						case 'Siemens'
							% Repeat the PLD for each TE first, then for repetitions:
							% Example 1: 4PLDs, 3 TEs, 2 repetitions -> 24 volumes [1 1 1 2 2 2 3 3 3 4 4 4 1 1 1 2 2 2 3 3 3 4 4 4]
							% Example 2: 4PLDs,  1 TE, 2 repetitions ->  8 volumes [1 2 3 4 1 2 3 4]
							jsonOut.PostLabelingDelay = repmat(jsonOut.PostLabelingDelay(:)', NumberEchoTimes, nVolumes/length(jsonOut.PostLabelingDelay)/NumberEchoTimes);
							jsonOut.PostLabelingDelay = jsonOut.PostLabelingDelay(:);
						case 'Philips'
							jsonOut.PostLabelingDelay = repmat(jsonOut.PostLabelingDelay(:)', 1, nVolumes/length(jsonOut.PostLabelingDelay));
							jsonOut.PostLabelingDelay = jsonOut.PostLabelingDelay(:);
					end
				else
					warning('Number of PLDs %d and TEs %d does not match the number of volumes %d\n', length(jsonOut.PostLabelingDelay), NumberEchoTimes, nVolumes);
				end
                
            end
        end
    end

    % Check for FME Hadamard sequences
    if isfield(jsonOut,'SeriesDescription')
        bTimeEncodedFME = ~isempty(regexp(jsonOut.SeriesDescription,'(Encoded_Images_Had)\d\d(_)\d\d(_TIs_)\d\d(_TEs)', 'once'));
        if bTimeEncodedFME
            startDetails = regexp(jsonOut.SeriesDescription,'\d\d(_)\d\d(_TIs_)\d\d(_TEs)', 'once');
            jsonOut.TimeEncodedMatrixSize = xASL_str2num(jsonOut.SeriesDescription(startDetails:startDetails+1));
            TimeEncodedInversionTimes = xASL_str2num(jsonOut.SeriesDescription(startDetails+3:startDetails+4));
            NumberEchoTimesFME = xASL_str2num(jsonOut.SeriesDescription(startDetails+10:startDetails+11));
			if NumberEchoTimesFME ~= NumberEchoTimes
				warning('TE count differs in DICOM %d and FME header %d\n', NumberEchoTimes, NumberEchoTimesFME);
			end
            fprintf('FME sequence, Hadamard-%d encoded images, %d TIs, %d TEs\n', jsonOut.TimeEncodedMatrixSize, TimeEncodedInversionTimes, NumberEchoTimes);
        end
    end
    
    % Check total acquired pairs for time encoded sequences
    if bTimeEncoded
		% At this stage, jsonOut.PostLabelingDelay has the repeated PLDs already.
		numberPLDs = jsonOut.TimeEncodedMatrixSize;
		% Determine the TotalAcquiredPairs
		jsonOut.TotalAcquiredPairs = nVolumes / (NumberEchoTimes * numberPLDs);
    end

end




