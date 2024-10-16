function jsonOut = xASL_bids_BIDSifyM0(jsonIn, jsonInASL, studyPar, pathM0In, pathM0Out, headerASL)
%xASL_bids_BIDSifyM0 Goes through the JSON and NII of an M0 file and makes sure that all the necessary conversions and checks to BIDS format are applied
%
% FORMAT: [jsonOutM0, jsonOutASL] = xASL_bids_BIDSifyM0(jsonIn, jsonASL, studyPar, pathM0In, pathM0Out, headerASL)
%
% INPUT:
%   jsonIn        - JSON with the input fields in BIDS format for M0 (REQUIRED)
%   jsonInASL     - JSON with the input fields in BIDS format for ASL (REQUIRED)
%   studyPar      - Manually entered study parameters (REQUIRED)
%   pathM0In      - Path to the source M0 file without extension (REQUIRED)
%   pathM0Out     - Path to the output M0 file without extension (REQUIRED)
%   headerASL     - NIfTI header of the ASL file (REQUIRED)
%
% OUTPUT: 
%   jsonOut       - Output JSON for M0
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    It makes all the conversions to a proper BIDS structure, checks the existence of all BIDS 
%                 fields, removes superfluous fields, checks all the conditions and orders the structure on 
%                 the output. It works according to the normal BIDS, or ASL-BIDS definition. It modifies 
%                 the NIfTI file to take into account several BIDS specifics.
%                 Specifically, it applies the previously calculated scalings.
%
% 1. Check the scaling in DICOMs
% 2. Check the JSON parameters
% 3. Save or move the NII to the correct location
%
% EXAMPLE: n/a
%
% __________________________________
% Copyright 2015-2023 ExploreASL

% Check if required fields exist in studyPar but not in jsonIn
jsonIn = xASL_bids_MergeStudyPar(jsonIn,studyPar,'m0');

% Create default output
jsonOut = jsonIn;

%% 1. Check the scaling in DICOMs
headerM0 = xASL_io_ReadNifti([pathM0In '.nii']);

if ~isempty(regexpi(jsonInASL.Manufacturer, 'Philips'))
    jsonOut.scaleFactor = xASL_adm_GetPhilipsScaling(jsonOut, headerM0);
else
    jsonOut.scaleFactor = 1;
end

%% 2. Check the JSON parameters
% Check echo time, for vectors
if isfield(jsonOut, 'EchoTime') && length(jsonOut.EchoTime)>1
    % Remove zero entries
    jsonOut.EchoTime = jsonOut.EchoTime(jsonOut.EchoTime ~= 0);
end

if isfield(studyPar, 'TotalReadoutTime')
    jsonOut.TotalReadoutTime = studyPar.TotalReadoutTime;
end

if isfield(jsonInASL, 'SliceTiming')
    % Issue a warning if the SliceTiming was already existing for M0, but still overwrite with ASL one
    if isfield(jsonOut, 'SliceTiming')
        warning('SliceTiming already existed for M0, overwriting with ASL');
    end
    
    if headerASL.dat.dim(3) == headerM0.dat.dim(3)
        % Either copy if the save number of slices in M0 as in ASL
        jsonOut.SliceTiming = jsonInASL.SliceTiming;
    else
        % Or recalculate for M0 if the number of slices differ
        jsonOut.SliceTiming = ((0:(headerM0.dat.dim(3)-1))')*(jsonInASL.SliceTiming(2)-jsonInASL.SliceTiming(1));
    end
else
    if isfield(jsonOut, 'SliceTiming')
        jsonOut = rmfield(jsonOut, 'SliceTiming');
        warning('Removing pre-existing SliceTiming from M0, as there was no SliceTiming for ASL');
    end
end

jsonOut.RepetitionTimePreparation = jsonOut.RepetitionTime;

% Check for Philips Look-Locker
% Look-Locker acquisition for Philips has often TR given as the length of a single readout of single volume, not as the entire cycle of all PLDs within a single labeling cycle
% The correct TR thus has to be correctly recalculated

% Verify that we are dealing with Philips Look-Locker sequence with flip angle below 90
if isfield(studyPar, 'LookLocker') && studyPar.LookLocker
	if isfield(jsonOut, 'Manufacturer') && strcmpi(jsonOut.Manufacturer, 'Philips')
		% Verify if TR is suspiciously low
		if isfield(jsonOut, 'FlipAngle') && jsonOut.FlipAngle < 90
			% If the flip angle is below 90 then also M0 was acquired as Look-Locker and would typically be as long as the PLD
			if isfield(studyPar, 'PostLabelingDelay') && isfield(jsonOut, 'RepetitionTimePreparation') && max(studyPar.PostLabelingDelay) > max(jsonOut.RepetitionTimePreparation)
				numberVolumes = size(headerM0.dat, 4);
				if isfield(jsonOut, 'PhilipsNumberTemporalScans') && max(jsonOut.PhilipsNumberTemporalScans) > 1
					% PhilipsNumberTemporalScans gives out the number of scans, so if it is divisible by the total number of volumes, then we can create the TR vector
					numberPhases = max(jsonOut.PhilipsNumberTemporalScans);
				elseif isfield(studyPar, 'ArterialSpinLabelingType') && strcmpi(studyPar.ArterialSpinLabelingType, 'PCASL') && isfield(studyPar, 'LabelingDuration')
					% If everything is provided in a single volume then the number of temporal scans is 1 and the number of volumes per cycle has to be quessed from the maximal PLD and labeling duration
					numberPhases = round(numberVolumes/(length(unique(studyPar.PostLabelingDelay)) + max(studyPar.LabelingDuration)/jsonOut.RepetitionTimePreparation));
					numberPhases = numberVolumes/numberPhases;
				elseif isfield(studyPar, 'LabelingType') && strcmpi(studyPar.LabelingType, 'PCASL') && isfield(studyPar, 'LabelingDuration')
					numberPhases = round(numberVolumes/(length(unique(studyPar.PostLabelingDelay)) + max(studyPar.LabelingDuration)/jsonOut.RepetitionTimePreparation));
					numberPhases = numberVolumes/numberPhases;
				else
					% Or maximum PLD only for PASL
					numberPhases = round(numberVolumes/(length(unique(studyPar.PostLabelingDelay))));
					numberPhases = numberVolumes/numberPhases;
				end
				% Create the TR vector
				if numberPhases && mod(numberVolumes, numberPhases) == 0
					jsonOut.RepetitionTimePreparation = (1:numberPhases) * jsonOut.RepetitionTimePreparation;
					jsonOut.RepetitionTimePreparation = repmat(jsonOut.RepetitionTimePreparation, [1 numberVolumes/numberPhases]);
				end
			end
		end
	end
end

%% 3. Save or move the NII to the correct location
xASL_adm_CreateDir(fileparts(pathM0Out));

% Validate the M0 output filename

[~,outputFileM0,outputExtensionM0] = xASL_fileparts([pathM0Out '.nii.gz']);
outputFilenameM0 = [outputFileM0 outputExtensionM0];
xASL_bids_ValidateNiftiName(outputFilenameM0,'m0scan');

% The NIfTI needs to be read and saved again
if (jsonOut.scaleFactor && jsonOut.scaleFactor~=1) || length(headerM0.dat.dim) < 4 || headerM0.dat.dim(4) == 1
    % Read NIfTI image
    imM0   = xASL_io_Nifti2Im([pathM0In '.nii']);
    
    % Apply the scaling
    if (jsonOut.scaleFactor && jsonOut.scaleFactor~=1)
        imM0 = imM0 .* jsonOut.scaleFactor;
    end
    
    % Save the NIfTI to a new location
    xASL_io_SaveNifti([pathM0In '.nii'],[pathM0Out '.nii.gz'],imM0,[],1,[]);
    
    % Delete original Nifti if not the same file
    if ~strcmp(pathM0In, pathM0Out)
        xASL_delete([pathM0In '.nii']);
    end
else
    % Move the M0
    xASL_Move([pathM0In '.nii'], [pathM0Out '.nii.gz'],1);
end

jsonOut = rmfield(jsonOut,'scaleFactor');


end