function [CorrectedImage, x] = xASL_im_LookLockerIntensityCorrection(image, x)
%xASL_im_LookLockerIntensityCorrection Correct Look-Locker 4D image for steady state magnetisation effects
%
% FORMAT: x = xASL_im_LookLockerIntensityCorrection(x)
%
% INPUT:
%   4Dimage     - Look-Locker 4D image to be corrected.
%   x           - struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%   ImageSavePath  - Path of corrected image to be saved to
% OUTPUT:
%   x           - struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%                 assigns the following values:
%   Corrected4Dimage  - Look-Locker 4D image corrected for steady state magnetisation effects
%-----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function corrects Look-Locker PWI4D for steady state magnetisation effects

% check if PWI or other 4D image
NUniquePLD = numel(unique(x.Q.Initial_PLD));

if size(image,4) > NUniquePLD % 4D non-PWI image, assuming control for Look-Locker M0 creation
    disp('Correcting Look-Locker control images for steady state magnetisation effects')
    for iAverage = 1 : (size(image,4)/NUniquePLD) % number of PLDs
        for iPLD = 1 : NUniquePLD % number of PLDs
            CorrectedImage(:,:,:,(iPLD +(NUniquePLD * (iAverage-1)))) = image(:,:,:,(iPLD + (NUniquePLD * (iAverage-1))))/((cos(2*pi/360*x.Q.FlipAngle))^(iPLD-1)); % correct control images for Look-Locker readout per PLD number, per average, see Günther, M., Bock, M. and Schad, L.R. (2001), Arterial spin labeling in combination with a look-locker sampling strategy: Inflow turbo-sampling EPI-FAIR (ITS-FAIR). Magn. Reson. Med., 46: 974-984. https://doi.org/10.1002/mrm.1284
        end
    end
else % 4D PWI images
    disp('Correcting Look-Locker PWI images for steady state magnetisation effects')
    for iPLD = 1 : size(unique(x.Q.Initial_PLD),1) % number of PLDs
        CorrectedImage(:,:,:,iPLD) = image(:,:,:,iPLD)/((cos(2*pi/360*x.Q.FlipAngle))^(iPLD-1)); % correct PWI for Look-Locker readout per PLD number, see Günther, M., Bock, M. and Schad, L.R. (2001), Arterial spin labeling in combination with a look-locker sampling strategy: Inflow turbo-sampling EPI-FAIR (ITS-FAIR). Magn. Reson. Med., 46: 974-984. https://doi.org/10.1002/mrm.1284
    end
end
end