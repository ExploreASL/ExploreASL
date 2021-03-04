function SliceTime = xASL_quant_SliceTimeVector(x,imPathMatrix)
%xASL_quant_SliceTimeVector Calculates the vector of slice-times for the given number of slices
%
% FORMAT: SliceTime = xASL_quant_SliceTimeVector(x,imPathMatrix)
%
% INPUT:
%   x            - struct containing pipeline environment parameters (REQUIRED)
%   imPathMatrix - path to the image or an image matrix used to calculated the number of slices to calculate the SliceTime for (REQUIRED)
% OUTPUT:
%   SliceTime  - a vector of SliceTimes with the length corresponding to the third dimension in imPathMatrix
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function takes the x.SliceReadoutTime, which can be a vector or scalar and creates a vector
%              out of it with the correct length. It also checks the x.readout_dim for which it returns 0.
%              It loads the image from imPathMatrix and calculates the SliceTime according to the number of slices in the third dimension
%              If a path is given, it also checks if it can find a JSON and loads it a looks for SliceTime as a replacement for the value in x-struct
%
% 0. Admin
% 1. ShortestTR
% 2. Assign the vector value
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: SliceTime = xASL_quant_SliceTimeVector(x,32)
% __________________________________
% Copyright 2015-2021 ExploreASL

%% -----------------------------------------------------------------------------------------------------------------------------------------------------
%% 0.Admin
if nargin < 1 || isempty(x)
	error('The x-structure needs to be provided');
end

if nargin < 2 || isempty(imPathMatrix) 
	error('The imPathMatrix parameter needs to be provided');
end

% Loads the matrix if path is provided
if ischar(imPathMatrix)
	% If path was given, load and check its third dimension
	imMatrix = xASL_io_ReadNifti(imPathMatrix);
	nSlices = size(imMatrix.dat,3);
	
	% Additionally, look for the JSON sidecar and its content
	[Fpath, Ffile, ~] = xASL_fileparts(imPathMatrix);
	pathJson = fullfile(Fpath,[Ffile,'.json']);
	
	% There might be an additional 'r' at the start to remove
	if ~xASL_exist(pathJson,'file') && (Ffile(1) == 'r')
		pathJson = fullfile(Fpath,[Ffile(2:end),'.json']);
	end
	
	% Check if the JSON exists
	if xASL_exist(pathJson,'file')
		parmsLocal = spm_jsonread(pathJson);
		% Loads and converts to the Legacy format
		parmsLocal = xASL_bids_parms2BIDS([], parmsLocal, 0);
		
		% Gets the SliceReadoutTime out of it
		if isfield(parmsLocal,'Q') && isfield(parmsLocal.Q,'SliceReadoutTime') && ~isempty(parmsLocal.Q.SliceReadoutTime)
			% If the field from JSON differs with the x-struct, but the JSON variant has the correct length, then replace it
			if ~isequal(parmsLocal.Q.SliceReadoutTime(:),x.Q.SliceReadoutTime(:)) &&...
					( (length(parmsLocal.Q.SliceReadoutTime) == 1) || (length(parmsLocal.Q.SliceReadoutTime) == nSlices))
				x.Q.SliceReadoutTime = parmsLocal.Q.SliceReadoutTime;
				warning(['Replacing the x-struct SliceReadoutTime with the one from: ' pathJson]);
			end
		end
	end
else
	% If a matrix was given, check its third dimension
	nSlices = size(imPathMatrix,3);
end

% The readout_dim needs to be provided
if ~isfield(x,'readout_dim')
	error('x.readout_dim field is missing');
end

if strcmpi(x.readout_dim, '3D')
	% For 3D sequences, zero is returned
	SliceTime = 0;
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
x = xASL_quant_SliceReadoutTime_ShortestTR(x);

%% -----------------------------------------------------------------------------------------------------------------------------------------------------
%% 2. Assign the vector value
% For 2D it either uses the vector or it replicates the scalar to the correct length
% Non-2D cases were solved in the admin

if length(x.Q.SliceReadoutTime) == 1
	SliceTime = (0:1:(nSlices-1)) * x.Q.SliceReadoutTime;
elseif length(x.Q.SliceReadoutTime) == nSlices
	SliceTime = x.Q.SliceReadoutTime;
else
	error('x.Q.SliceReadoutTime has to be a scalar or match the number of slices');
end

end
