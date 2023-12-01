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
if x.modules.asl.bTimeEncoded
    % Here we select the 1st TE and the control images
    % We thus calculate the size of each Hadamard block as the number of Hadamard phases and TEs
    blockSize = x.Q.TimeEncodedMatrixSize * x.Q.NumberEchoTimes;
    
    % For example for 64 volumes and 2 repetitions with 8 PLDs and 4 TEs, it takes volume 1 and 33
    imMeanControl = imASLTimeSeries(:,:,:,1:blockSize:end);
    
elseif x.modules.asl.bMultiPLD
    % Create mean control in native space
    if x.modules.asl.bContainsDeltaM
        imMeanControl = imASLTimeSeries;
        Initial_PLD = x.Q.Initial_PLD;
    elseif (isfield(x.Q,'LookLocker') && x.Q.LookLocker) || x.modules.asl.bContainsDeltaM
        imMeanControl = xASL_quant_GetControlLabelOrder(imASLTimeSeries);
        % Get the PLD for the given sequence by dividing with two
        Initial_PLD = x.Q.Initial_PLD(1:(length(x.Q.Initial_PLD)/2));
    else
        imMeanControl = xASL_quant_GetControlLabelOrder(imASLTimeSeries);
        % Get the PLD for the given sequence
        Initial_PLD = x.Q.Initial_PLD(1:2:end);
    end
    
    % Get unique PLDs
    idealPLD = unique(Initial_PLD);
    
    % Find the index of the one closest to 2000 ms
    [~, ind] = min(abs(idealPLD-2000));
    
    % Pick up the ideal PLD as the one closest to 2000 ms
    idealPLD = idealPLD(ind(1));
    
    if (isfield(x.Q,'LookLocker') && x.Q.LookLocker) || x.modules.asl.bContainsDeltaM
        % For Look-Locker, get the last one
        NUniquePLD = numel(unique(Initial_PLD));
        idealPLDLookLocker = unique(Initial_PLD);
        idealPLD = idealPLDLookLocker(end);
    else
        % For normal mutli-PLD, get the latest PLD
        idealPLD = idealPLD(end);
    end
    
    % Find all dynamics with that PLD
    if  isfield(x.Q,'LookLocker') && x.Q.LookLocker
        disp('Correcting Look-Locker control images for steady state magnetisation effects')
        for iAverage = 1 : (size(imMeanControl,4)/ NUniquePLD) % number of PLDs
            for iPLD = 1 :  NUniquePLD % number of PLDs
                imMeanControlCorrected(:,:,:,(iPLD +(NUniquePLD * (iAverage-1)))) = imMeanControl(:,:,:,(iPLD + (NUniquePLD * (iAverage-1))))/((cos(2*pi/360*x.Q.FlipAngle))^(iPLD-1)); % correct control images for Look-Locker readout per PLD number, per average, see Günther, M., Bock, M. and Schad, L.R. (2001), Arterial spin labeling in combination with a look-locker sampling strategy: Inflow turbo-sampling EPI-FAIR (ITS-FAIR). Magn. Reson. Med., 46: 974-984. https://doi.org/10.1002/mrm.1284
            end
        end
    end
    imMeanControl = imMeanControlCorrected(:,:,:,Initial_PLD==idealPLD);
    
else
    % Create mean control in native space
    imMeanControl = xASL_quant_GetControlLabelOrder(imASLTimeSeries);
    
end

%% Calculate mean over the repetitions
imMeanControl = xASL_stat_MeanNan(imMeanControl, 4);

end
