function xASL_im_CleanupWMHnoise(InputPath, OutputPath, MinLesionVolume, pThresh)
%xASL_im_CleanupWMHnoise Threshold white matter lesions,
% acknowledging the fact that they may be confluent with subresolution connection
% through a dilation. This part is executed conservatively, as FLAIR hyperintensities
% inside the GM can be erroneously segmented as WMH, and should not be lesion-filled
% (otherwise these cannot be fixed later in the Structural module)
%
% Note that LST lesion filling expects a probability map, doesnt work nicely with binary mask

%% ---------------------------------------------------------------------
%% Admin

if ~exist('MinLesionVolume', 'var')
    MinLesionVolume     = 15; % default
end
if ~exist('pThresh', 'var')
    pThresh             = 0.5; % default
end

% Determine minimum volumes
SphereVolConst          = (4/3)*pi;
VolumeR                 = (MinLesionVolume/SphereVolConst)^(1/3)+0.50; % after dilation with 1 (note R = single side)
MinVolumeDilated        = VolumeR^3*SphereVolConst;



%% ---------------------------------------------------------------------
%% Load NIfTI
volWMH                  = xASL_io_ReadNifti(InputPath);
imWMH                   = xASL_io_Nifti2Im(InputPath);

% Calculate voxel size and volume
voxelVolumeWMH          = prod(sqrt(sum(volWMH.mat(1:3,1:3).^2)));

fprintf('%s','Cleaning up WMH_SEGM:  0%');



%% ---------------------------------------------------------------------
%% Run the thresholding

% Threshold at low value not to miss any connections
imWMHthresholded        = imWMH > pThresh;
OutputImage = imWMH;

% First dilate
imThreshDilated         = xASL_im_DilateErodeFull( imWMHthresholded, 'dilate', xASL_im_DilateErodeSphere(2));

% Label the connected components
[imLabel, numLabel]     = spm_bwlabel(double(imThreshDilated),26);

% Go through all the components of the lesion image
for ii = 1:numLabel
    xASL_TrackProgress(ii, numLabel);
    tL                  = imWMH(imLabel==ii);
    LesionVolume        = voxelVolumeWMH * xASL_stat_SumNan(tL(:));

    % Remove volumes that are too small
    if  LesionVolume < MinVolumeDilated
        imWMHthresholded(imLabel == ii) = 0;
        OutputImage(imLabel == ii) = 0; 
        % LST lesion filling expects a probability map, doesnt work nicely with binary mask
    end
end

fprintf('\n');

xASL_io_SaveNifti(InputPath, OutputPath, OutputImage, [], false);

end
