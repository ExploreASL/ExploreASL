function SliceTime = xASL_quant_SliceTimeVector(x,nSlices)
%xASL_quant_SliceTimeVector Calculates the vector of slice-times for the given number of slices
%
% FORMAT: SliceTime = xASL_quant_SliceTimeVector(x,nSlices)
%
% INPUT:
%   x          - struct containing pipeline environment parameters (REQUIRED)
%   nSlices    - number of slices in ASL native space to calculate the SliceTime for (REQUIRED)
% OUTPUT:
%   SliceTime  - a vector of SliceTimes with the length nSlices
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function takes the x.SliceReadoutTime, which can be a vector or scalar and creates a vector
%              out of it with the correct length. It also checks the x.readout_dim for which it returns 0.
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

if nargin < 2 || isempty(nSlices) 
	error('The nSlices parameter needs to be provided');
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
