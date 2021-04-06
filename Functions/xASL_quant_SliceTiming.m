function SliceTiming = xASL_quant_SliceTiming(x, inputIm)
%xASL_quant_SliceTiming Takes the SliceReadoutTime (xASL-legacy or BIDS) and calculates a vector with the timings of slice readouts relative 
% to the start of readout (BIDS definition)
%
% FORMAT: SliceTiming = xASL_quant_SliceTiming(x, inputIm)
%
% INPUT:
%   x            - struct containing pipeline environment parameters (REQUIRED)
%   inputIm      - path to the image or an image matrix used to calculated the number of slices to calculate the SliceTiming for (REQUIRED)
% OUTPUT:
%   SliceTiming  - a vector of slice times defining the start of readout per slices with respect to the start of the readout of the first slices
%                  with the length corresponding to the third dimension in inputIm
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function takes the x.Q.SliceReadoutTime, which can be:
%              i) a vector of slice readout-timings relative to the start of readout of the first slice or (i.e. the BIDS definition of SliceTiming)
%              ii) a scalar with difference in readout times between the consecutives slices (i.e. the xASL legacy definition of SliceTiming)
%              The function creates a vector (of the relatives timings for each slices) out of it with the correct length corresponding to the number of slices in the inputIm 
%              corresponding to the BIDS definition. It also checks the x.readout_dim, and for 3D readouts it returns 0.
%              It loads the image from inputIm and calculates the SliceTiming according to the number of slices in the third dimension
%              If a path is given, it also checks if it can find a JSON sidecar, then it loads the JSON sidecar, and looks for SliceTiming inside it. If
%              SliceTiming/SliceReadoutTime is found in the JSON sidecar, it prioritize it over the value in the x-struct
%
% 0. Admin
% 1. ShortestTR
% 2. Assign the vector value and check for vector consistency
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE 1: SliceTiming = xASL_quant_SliceTiming(x, 'ASL4D.nii')
%
% EXAMPLE 2: x.Q.SliceReadoutTime = [0 30 60 90 120]
%            SliceTiming = xASL_quant_SliceTiming(x, 'ASL4D.nii')
%
% EXAMPLE 3: x.Q.SliceReadoutTime = [30]
%            SliceTiming = xASL_quant_SliceTiming(x, 'ASL4D.nii')
% __________________________________
% Copyright 2015-2021 ExploreASL

%% -----------------------------------------------------------------------------------------------------------------------------------------------------
%% 0.Admin
if nargin < 1 || isempty(x)
	error('The x-structure needs to be provided');
end

if nargin < 2 || isempty(inputIm) 
	error('The inputIm parameter needs to be provided');
end

% If a path is provided, then it needs to load the image to obtain its number of dimensions and it can also check the file JSON sidecar
% to see if a different SliceReadoutTime is no provided in the JSON than in the x-struct
if ischar(inputIm)
	% Loads the image information and obtains the number of slices
	imMatrix = xASL_io_ReadNifti(inputIm);
	nSlices = size(imMatrix.dat, 3);
	
	% Additionally, look for the JSON sidecar of the image file
	[Fpath, Ffile, ~] = xASL_fileparts(inputIm);
	pathJson = fullfile(Fpath, [Ffile, '.json']);
	
	% We need to look for the sidecar of the original file in case a resampled file was provided (with a prefix 'r')
	if ~xASL_exist(pathJson, 'file') && (Ffile(1) == 'r')
		pathJson = fullfile(Fpath, [Ffile(2:end), '.json']);
	end
	
	% Check if the JSON sidecar of the original file exists
	if xASL_exist(pathJson, 'file')
		parmsLocal = spm_jsonread(pathJson);
		% Loads the JSON file and converts it to the Legacy format
		parmsLocal = xASL_bids_parms2BIDS([], parmsLocal, 0);
		
		% Gets the SliceReadoutTime out of it
		if isfield(parmsLocal, 'Q') && isfield(parmsLocal.Q, 'SliceReadoutTime') && ~isempty(parmsLocal.Q.SliceReadoutTime)
			% If SliceReadoutTime differs from the x-struct version, but the JSON variant has the correct length, then replace it
			if ~isequal(parmsLocal.Q.SliceReadoutTime(:), x.Q.SliceReadoutTime(:)) &&...
					( (length(parmsLocal.Q.SliceReadoutTime) == 1) || (length(parmsLocal.Q.SliceReadoutTime) == nSlices))
				x.Q.SliceReadoutTime = parmsLocal.Q.SliceReadoutTime;
				warning(['Replacing the x-struct SliceReadoutTime with the one from: ' pathJson]);
			end
		end
	end
else
	% If an image matrix is given, the directly check the number of slices
	nSlices = size(inputIm, 3);
end

% The readout_dim needs to be provided
if ~isfield(x, 'readout_dim')
	error('x.readout_dim field is missing');
end

if strcmpi(x.readout_dim, '3D')
	% For 3D sequences, zero is returned
	SliceTiming = 0;
	return;
end

if ~strcmpi(x.readout_dim, '2D') 
	% It only works with 2D and 3D
	error(['Unknown x.readout_dim value:' x.readout_dim]);
end

% The SliceReadoutTime needs to be provided
if ~isfield(x, 'Q')
	error('x.Q field missing');
elseif ~isfield(x.Q, 'SliceReadoutTime') || isempty(x.Q.SliceReadoutTime)
	error('x.Q.SliceReadoutTime missing or invalid ');
end

%% -----------------------------------------------------------------------------------------------------------------------------------------------------
%% 1. ShortestTR
% If SliceReadoutTiem is specified as "shortestTR", it calculates it with the knowledge of TR and PLD
% If a scalar or vector is given for SliceReadoutTime, then this function doesn't do anything
x = xASL_quant_SliceTiming_ShortestTR(x);

%% -----------------------------------------------------------------------------------------------------------------------------------------------------
%% 2. Assign the vector value
% For 2D it either uses the vector or it replicates the scalar to the correct length
% Non-2D cases were solved in the admin

if length(x.Q.SliceReadoutTime) == 1
	SliceTiming = (0:1:(nSlices-1)) * x.Q.SliceReadoutTime;
elseif length(x.Q.SliceReadoutTime) == nSlices
	SliceTiming = x.Q.SliceReadoutTime;
	
	% Check for SliceTiming consistency between slices if a vector was provided
	if nSlices > 1
		% SliceTimingDiff is a temporary variable that describes the difference between neighboring slices, 
		% has the same length as the number of slices - 1. And can be negative for Multi-Band. It is used only internally
		% to determine if the SliceTiming is consistent and/or multi-band
		SliceTimingDiff = SliceTiming(2:end) - SliceTiming(1:(end-1));
		
		% See if the slice-timing difference between slices is equal for all slices
		[SliceTimingDiffUnique,SliceTimingDiffIndex] = uniquetol(SliceTimingDiff,0.01);
		
		% Make sure that the positive values are listed first
		[SliceTimingDiffUnique,SliceTimingDiffIndexSort] = sort(SliceTimingDiffUnique,'descend');
		SliceTimingDiffIndex = SliceTimingDiffIndex(SliceTimingDiffIndexSort);
		
		if length(SliceTimingDiffUnique) > 2
			% If more than two values of between slice-timing difference are present, then the SliceTiming is potentially inconsistent
			fprintf('Warning: Inconsistent SliceTiming between slices\n');
		elseif length(SliceTimingDiffUnique) == 2
			% Two different values point to Multi-Band readout:
			% SliceTiming rises equally between slices (difference #1) and decreases between Multi-Band blocks (difference #2)
			bIsMultiBand = false;
			if SliceTimingDiffUnique(1) > 0 && SliceTimingDiffUnique(2) < 0
				% Multi-band should have a distinct pattern in SliceTiming [0 X 2*X 3*X .. 0 X 2*X 3*X .. 0]
				MultiBandFactor = SliceTimingDiffIndex(2);
				if mod(length(SliceTiming),MultiBandFactor) == 0
					% Try to reconstruct the pattern of MultiBand and check if that is really the case
					SliceTimingPattern = (0:(MultiBandFactor-1))*SliceTimingDiffUnique(1);
					SliceTimingPattern = repmat(SliceTimingPattern,[1 length(SliceTiming)/MultiBandFactor])';
					if ~sum(isnear(SliceTiming',SliceTimingPattern',5)==0)
						bIsMultiBand = true;
					end
				end
			end
			if bIsMultiBand
				fprintf('Warning: Multi-band pattern detected in SliceTiming\n');
			else
				fprintf('Warning: Inconsistent SliceTiming between slices\n');
			end
		end
	end
	
	
else
	error('x.Q.SliceReadoutTime has to be a scalar or match the number of slices');
end

end
