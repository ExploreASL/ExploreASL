function imMeanControl = xASL_quant_GetMeanControl(x, imASLTimeSeries)
%xASL_quant_GetControlLabelOrder Get control-label order of ASL time-series
%
% FORMAT: imMeanControl = xASL_quant_GetMeanControl(x, imASLTimeSeries)
%
% INPUT:
%   x               - structure containing fields with all information required to run this submodule (REQUIRED)
%   imASLTimeSeries - 4D matrix with multiple 3D ASL image volumes over time (REQUIRED)
%
% OUTPUT:
%   imMeanControl   - 3D ASL mean control image volume
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function calculates the mean control image. It automatically checks (and corrects if required)
%              the control order of ASL timeseries and takes into account if multiPLD or Hadamard data are provided
%
% EXAMPLE:
% imMeanControl = xASL_quant_GetMeanControl(x, imASLTimeSeries)
% 
% __________________________________
% Copyright (c) 2015-2021 ExploreASL

%% Admin
if nargin < 1 || isempty(x)
	error('Missing x input');
end

if nargin < 2 || isempty(imASLTimeSeries)
	error('Missing imASLTimeSeries input');
end

%% Select the control images based on the acquisition type
if isfield(x.modules.asl,'bTimeEncoded') && x.modules.asl.bTimeEncoded 
	% Here we select the 1st TE and the control images
	% We thus calculate the size of each Hadamard block as the number of Hadamard phases and TEs
	blockSize = x.Q.TimeEncodedMatrixSize * x.Q.NumberEchoTimes;
        
    % For example for 64 volumes and 2 repetitions with 8 PLDs and 4 TEs, it takes volume 1 and 33
    imMeanControl = imASLTimeSeries(:,:,:,1:blockSize:end);
else
	% Create mean control in native space (works also for multi-PLD)
	imMeanControl = xASL_quant_GetControlLabelOrder(imASLTimeSeries);
	
	% For multi-PLD, we need to select one PLD
	if isfield(x.modules.asl,'bMultiPLD') && x.modules.asl.bMultiPLD
		% Get the PLD for the given sequence
		Initial_PLD = x.Q.Initial_PLD(1:2:end);
		
		% Sort the unique PLDs and pick the ideal PLD
		idealPLD = sort(unique(Initial_PLD),'ascending');
		
		if x.Q.LookLocker
			% For Look-Locker, get the middle one
			idealPLD = idealPLD(round(numel(idealPLD)/2));
		else
			% For normal mutli-PLD, get the latest PLD
			idealPLD = idealPLD(end);
		end
		
		% Find all dynamics with that PLD
		imMeanControl = imMeanControl(:,:,:,Initial_PLD==idealPLD);
	end
end

%% Calculate mean over the repetitions
imMeanControl = xASL_stat_MeanNan(imMeanControl, 4);
	
end
