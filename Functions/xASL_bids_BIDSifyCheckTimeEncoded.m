function [jsonOut,bTimeEncoded,bTimeEncodedFME] = xASL_bids_BIDSifyCheckTimeEncoded(jsonIn,jsonOut)
%xASL_bids_BIDSifyCheckTimeEncoded Check for time encoded sequence
%
% FORMAT: [jsonOut,bTimeEncoded,bTimeEncodedFME] = xASL_bids_BIDSifyCheckTimeEncoded(jsonIn,jsonOut)
%
% INPUT:
%   jsonIn    - JSON with the input fields - from DICOMs (REQUIRED)
%   jsonOut   - Output JSON in progress from the parent function xASL_bids_BIDSifyASLJSON (REQUIRED)
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

    % Check for time encoded sequences
    if isfield(jsonOut,'TimeEncodedMatrixSize') && ~isempty(jsonOut.TimeEncodedMatrixSize) || ... % Should be 4, 8 or 12
        isfield(jsonOut,'TimeEncodedMatrixType') % Natural or walsh
        bTimeEncoded = true; 
    else 
        bTimeEncoded = false; 
    end


    % Check for specific time encoded sequence of FME (Fraunhofer Mevis)
    if isfield(jsonIn,'SeriesDescription')
        bTimeEncodedFME = ~isempty(regexp(jsonIn.SeriesDescription,'(Encoded_Images_Had)\d\d(_)\d\d(_TIs_)\d\d(_TEs)', 'once'));
    else
        bTimeEncodedFME = false;
    end

    if bTimeEncodedFME
        bTimeEncoded = true;    
    end

    % Default
    TimeEncodedEchoTimes = 1;
    TimeEncodedInversionTimes = 1;
    
    if bTimeEncoded
        if isfield(jsonOut,'EchoTime') && isfield(jsonOut,'PostLabelingDelay')
            % From the import, the length of EchoTime should correspond to the number of volumes
            if length(jsonOut.EchoTime)~=length(jsonOut.PostLabelingDelay)
                % So here, we first make sure that each PLD is repeated for the whole block of echo-times
                TimeEncodedEchoTimes = length(uniquetol(jsonOut.EchoTime,0.001)); % Obtain the number of echo times
                repeatedPLDs = repmat(jsonOut.PostLabelingDelay(:),1,TimeEncodedEchoTimes)';
                repeatedPLDs = repeatedPLDs(:);

                if length(repeatedPLDs) > length(jsonOut.EchoTime) || mod(length(jsonOut.EchoTime),length(repeatedPLDs)) ~= 0
                    warning('Did not succeed in repeating PLDs for each TE for Hadamard sequence import...');
                end
                % Make sure that number of volumes can be divided by the repeated PLDs
                jsonOut.PostLabelingDelay = repeatedPLDs;
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
            TimeEncodedEchoTimes = xASL_str2num(jsonOut.SeriesDescription(startDetails+10:startDetails+11));
            fprintf('FME sequence, Hadamard-%d encoded images, %d TIs, %d TEs\n', jsonOut.TimeEncodedMatrixSize, TimeEncodedInversionTimes, TimeEncodedEchoTimes);
            jsonOut.TimeEncodedEchoTimes = TimeEncodedEchoTimes;
        end
    end
    
    % Check total acquired pairs for time encoded sequences
    if bTimeEncoded
        if (TimeEncodedEchoTimes>1) && (TimeEncodedInversionTimes>1)
            % At this stage, jsonOut.PostLabelingDelay has the repeated PLDs already.
            numberPLDs = TimeEncodedInversionTimes + 1;
            % Determine the TotalAcquiredPairs
            NumberRepetitions = jsonOut.AcquisitionMatrix / (TimeEncodedEchoTimes * numberPLDs);
            jsonOut.TotalAcquiredPairs = jsonOut.TimeEncodedMatrixSize * NumberRepetitions;
        end
    end

end




