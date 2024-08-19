function [IMout] = xASL_im_FillNaNs(InputPath, UseMethod, bQuality, VoxelSize, bZeroes)
%xASL_im_FillNaNs Fill NaNs in image
% FORMAT: xASL_im_FillNaNs(InputPath[, UseMethod, bQuality, VoxelSize, bZeroes])
%
% INPUT:
%   InputPath   - path to image, or image matrix (REQUIRED)
%                 NB: if this is a path, it is overwritten
%   UseMethod   - one of the following methods to fill the NaNs:
%                 1: extrapolate by smoothing (useful for any image, e.g. extrapolating M0 smooth biasfield) (DEFAULT)
%                 2: set to 0 (useful for e.g. DARTEL displacement fields, where 0 means no displacement)
%                 3: extrapolate deformation field linearly - this method is adapted to SPM deformation fields XxYxZx3 or XxYxZx1x3
%                 (OPTIONAL, DEFAULT = 1)
%   bQuality    - affects kernel size of the extrapolation performed in method 1
%                 (extrapolating by smoothing goes much faster with smaller
%                 kernel but this creates ringing artifacts)
%                 (OPTIONAL, DEFAULT = 1)
%   VoxelSize   - [X Y Z] vector for voxel size 
%                 (OPTIONAL, DEFAULT = [1.5 1.5 1.5])
%   bZeroes     - boolean for converting zeroes to NaNs and also changing them accordingly (OPTIONAL, DEFAULT = false)
%
%
% OUTPUT:
%   IMout       - image matrix in which NaNs were filled
% --------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function fills any NaNs in an image. In SPM, any voxels
%              outside the boundary box/field of view are filled by NaNs
%              when resampling. These NaNs can confuse some algorithms,
%              hence it doesn't hurt replacing them in some cases (e.g. for
%              flowfields). Also, smoothing restricted in a mask is done in
%              ExploreASL with the function xASL_im_ndnanfilter, after
%              first setting all voxels outside the mask to NaN. In this
%              case, this functon can be useful to extrapolate the smoothed
%              image to avoid any division artifact near brain edges (e.g.
%              for reducing the M0 image to a smooth biasfield).
%              This function performs the following 3 steps:
%
%              1. Load image
%              2. Replace NaNs using any of the methods
%               Method 1: extrapolate smooth over NaNs:
%                         Use this for cases where NaNs need to get border values (e.g. extrapolating a masked image)
%               Method 2: set NaNs to zeroes
%               Method 3: Use this for cases where NaNs need to expand linearly (e.g., extrapolating flowfields)
%              3. Save image
%
% --------------------------------------------------------------------------------------------------------------------
% EXAMPLE:   for filling NaNs in an image: xASL_im_FillNaNs('/MyStudy/sub-001/ASL_1/M0.nii');
% EXAMPLE2:  for fixing flowfield edges: xASL_im_FillNaNs('/MyStudy/sub-001/y_T1.nii', 3);
% __________________________________
% Copyright (C) 2015-2024 ExploreASL

%% Admin
IM = NaN; % default

if nargin<1 || isempty(InputPath)
    warning('Illegal input path, skipping');
    return;
end

if nargin<2 || isempty(UseMethod)
    UseMethod = 1;
elseif UseMethod<0 || UseMethod>3
    warning('Illegal option chosen, setting to default value method 1');
    UseMethod = 1;
end

if nargin<3 || isempty(bQuality)
    bQuality = 1; % default
elseif bQuality~=0 && bQuality~=1
    warning('Illegal bQuality setting!!!');
    bQuality = 1; % default
end

if UseMethod == 1 || UseMethod == 3 % only methods that needs VoxelSize
    if nargin < 4 || isempty(VoxelSize)
        if ~isnumeric(InputPath) && xASL_exist(InputPath, 'file')
            nii = xASL_io_ReadNifti(InputPath);
            VoxelSize = nii.hdr.pixdim(2:4);
		elseif UseMethod == 3
			% The Method 3 also need voxel-size, but it can estimate the voxel-size from the deformation fields
			% Therefore in case of missing VoxelSize, it is initialized as an empty array and this is automatically
			% estimated later when the image is properly loaded
			VoxelSize = [];
        elseif bQuality
			warning('Voxel size not specified, using default [1.5 1.5 1.5]');
			% only when high quality will we use this parameter,
			% so we dont give this warning for low quality
			VoxelSize = [1.5 1.5 1.5];
        else
			VoxelSize = [1.5 1.5 1.5];
        end
    end
else
    % We don't want the warning for methods other than 1, and we won't use VoxelSize for methods other than 1, but we predefine it here anyway.
    VoxelSize = [1.5 1.5 1.5];
end

if nargin<5 || isempty(bZeroes)
    if UseMethod == 3
        bZeroes = true;
    else
        bZeroes = false;
    end
end

%% ------------------------------------------------------------------------------
%% 1. Load image
IM = xASL_io_Nifti2Im(InputPath); % allows loading image from path or from image

if bZeroes
    % Zeros are also converted to NaNs - because this is sometimes unclear in the flow fields
    IM(IM==0) = NaN;
end

% Also catch Infs or imaginary numbers or whatever. Zeroes are already catched above in the main function.
IM(~isfinite(IM)) = NaN;

%% ------------------------------------------------------------------------------
%% 2. Replace NaNs using any of the methods
% Run only if there were actual NaNs present
if sum(isnan(IM(:))) > 0
	switch UseMethod
		case 1
			% Method 1 is run independently for each 3D volume
			fprintf('%s', 'Cleaning up NaNs in image:  0%');
			TotalDim = size(IM,4)*size(IM,5)*size(IM,6)*size(IM,7);
			CurrentIt = 1;
			for i4=1:size(IM,4)
				for i5=1:size(IM,5)
					for i6=1:size(IM,6)
						for i7=1:size(IM,7)
							xASL_TrackProgress(CurrentIt, TotalDim);

							% Use this for cases where NaNs need to get border values (e.g. extrapolating a masked image)
							IM(:,:,:,i4,i5,i6,i7) = xASL_im_ExtrapolateSmoothOverNaNs(IM(:,:,:,i4,i5,i6,i7), bQuality, VoxelSize);
							CurrentIt = CurrentIt+1;
						end
					end
				end
			end
			fprintf('\n');
			IMout = IM;
		case 2
			IM(isnan(IM)) = 0;
			IMout = IM;
		case 3
			% Use this for cases where NaNs need to expand linearly (e.g., extrapolating deformation flowfields)
			% First handle the two input options XxYxZx1x3 and XxYxZx3 by always changing it to XxYxZx3
			if ndims(IM) == 5 && size(IM, 4) == 1 && size(IM, 5) == 3
				% Squeeze from XxYxZx1x3 to XxYxZx3
				bAddDimension = true;
				IM = squeeze(IM);
			elseif ndims(IM) == 4 && size(IM, 4) == 3
				% No need to squeeze dimensions
				bAddDimension = false;
			else
				error('Selected method works only for 4D or 5D deformation fields of sizes XxYxZx3 or XxYxZx1x3');
			end
			IM = xASL_im_ExtrapolateLinearlyOverNaNs(IM, [], VoxelSize);

			% Recreate the correct dimensions for the output
			if bAddDimension
				IMout(:, :, :, 1, 1:3) = IM; 
			else
				IMout = IM;
			end
		otherwise
			error('Invalid option');
	end
else
	IMout = IM;
end

%% ------------------------------------------------------------------------------
%% 3. Save image
% Here we only output an image if the input path was a path and not an
% image matrix
if ~isnumeric(InputPath) && xASL_exist(InputPath, 'file')
    xASL_io_SaveNifti(InputPath, InputPath, IMout, [], 0);
end


end

%% ========================================================================================
%% ========================================================================================
%% METHOD 1
function [IM] = xASL_im_ExtrapolateSmoothOverNaNs(IM, bQuality, VoxelSize)
%xASL_im_ExtrapolateLinearlyOverNaNs

    if sum(isnan(IM(:))) < numel(IM)
        while sum(sum(sum(isnan(IM))))>0 % loop extrapolation until full FoV is filled
            if bQuality==1
                IM = xASL_im_ndnanfilter(IM,'gauss', double([12 12 12]./VoxelSize), 2); % 2 at the end means extrapolation only
            else
                IM = xASL_im_ndnanfilter(IM,'gauss', double([2 2 2]), 2); % smaller kernel goes much faster
            end
        end
    else
        % Image contains only NaNs - skip this image
	end
end

%% ========================================================================================
%% ========================================================================================
%% METHOD 3
function [IM] = xASL_im_ExtrapolateLinearlyOverNaNs(IM, nEdgeVoxels, VoxelSize)
%xASL_im_ExtrapolateLinearlyOverNaNs
% The input of this function is 3D (iterations of other dimensions are done
% in the main function above)
%
% IM          - input image
% nEdgeVoxels - number of voxels at the edges of the non-nan values that we check for weird interpolation values, 
% to be set to NaN and ignored in the extrapolation
% VoxelSize   - voxel size is needed for extrapolating correctly as the transformation is given in mm and not voxels

%% 0. Admin
if nargin<2 || isempty(nEdgeVoxels)
    nEdgeVoxels = 0;
    % We can try this later, but it seems to work less well than doing the same based on the distance transform
end

if nargin < 3 || isempty(VoxelSize)
	% The voxel-size is not provided and needs to be estimated from the transformation maps
	% We set it here to empty in admin and estimate it later
	VoxelSize = [];
end

%% Create masks and estimate voxel size
% Create NaN mask for all 3D volumes together
maskNaN = sum(isnan(IM), 4) > 0;

% The voxel-size needs to be estimated
if isempty(VoxelSize)
	% For each volume, we calculate the difference of neighboring voxels
	% We calculate the mean over nonNaN values, which gives the average update including the correct direction
	diffIM = IM(2:end, :, :, 1) - IM(1:end-1, :, :, 1);
	diffIM = diffIM(~maskNaN(2:end, :, :));
	VoxelSize(1) = xASL_stat_MeanNan(diffIM);


	diffIM = IM(:, 2:end, :, 2)-IM(:, 1:end-1, :, 2);
	diffIM = diffIM(~maskNaN(:, 2:end, :));
	VoxelSize(2) = xASL_stat_MeanNan(diffIM);

	diffIM = IM(:, :, 2:end, 3)-IM(:, :, 1:end-1, 3);
	diffIM = diffIM(~maskNaN(:, :, 2:end));
	VoxelSize(3) = xASL_stat_MeanNan(diffIM);
end

%% 0.5 Fix interpolation edges
% The first layer(s) next to the NaNs are sometimes interpolated as well
% So here we try to estimate which voxels those are, by looking at the distance from the NaNs
% together with the transformation difference from a smoothed transformation field.

% First we get the distance transform from mask (1 where NaNs are) inwards
distMapInward = xASL_im_DistanceTransform(maskNaN);

% Then we create a copy of the image that is smoothed (outside of NaNs)
IMsmooth = IM;
for iDim = 1:3
	IMsmooth(:, :, :, iDim) = xASL_im_ndnanfilter(IM(:, :, :, iDim), 'gauss', [8 8 8], 1);
end

% We obtain an absolute difference image between the original transformation field, and the smoothed transformation field
% This difference will show us how "extreme" the transformations are.
deltaIm = mean(abs(IM-IMsmooth), 4);
% This absolute difference image, is now divided by the distance.
% Because we want to detect extreme transformation values that are close to the border (low distance)
% So this would result in a high relativeDelta
relativeDelta = deltaIm./distMapInward;

% The relativeDelta is Rician-ish distributed, so we use the median and MAD to obtain a threshold for the relativeDelta
medianRelativeDelta = xASL_stat_MedianNan(relativeDelta(:));
MADRelativeDelta = xASL_stat_MadNan(relativeDelta(:));
thresholdIs = medianRelativeDelta + 3.5*MADRelativeDelta;

% We want to have a hard cutoff of distance from NaNs that we won't change
% This is 10% of the transformation field size (e.g, for 121*145*121 transformation field, this is 13 voxels)
sizeIm = size(IM);
if numel(sizeIm)<3
    warning('Something wrong with this flowfield');
end
averageSize = mean(sizeIm(1:3));
edgeDistance = round(0.1*averageSize);

% Set the voxels at the border where extreme values were detected to NaN as well
% We update the NaN mask instead of modifying the image
maskNaN(relativeDelta > thresholdIs & distMapInward < edgeDistance) = 1;

%% Extrapolation of values based on the edge values 
% Calculate the relative distance map from the non-Nan towards NaN. Note that for non-NaNs, these relative coordinates are zero
[~, distX, distY, distZ] = xASL_im_DistanceTransform(~maskNaN);

% We calculate the meshgrid - three 3D matrices that give, respectively, for each voxel its X, Y, and Z coordinates.
[gridX, gridY, gridZ] = ndgrid(1:size(IM, 1), 1:size(IM, 2), 1:size(IM, 3));

% We calculate the coordinates of the border voxel, note that for non-NaNs this gives the original voxel location
borderX = distX + gridX;
borderY = distY + gridY;
borderZ = distZ + gridZ;

for iDim = 1:3
	% To each NaN voxel, we assign the border value in the original image
	% Note that for non-NaN voxels, this keeps the exact original value
	IM(:, :, :, iDim) = IM(borderX + size(IM, 1) * (borderY-1 + size(IM, 2) * (borderZ-1 + size(IM, 3) * (iDim-1))));
end

% We now subtract the voxel-size adjusted relative coordinate difference from these newly filled voxels
IM(:, :, :, 1) = IM(:, :, :, 1) + distX(:, :, :)*VoxelSize(1);
IM(:, :, :, 2) = IM(:, :, :, 2) + distY(:, :, :)*VoxelSize(2);
IM(:, :, :, 3) = IM(:, :, :, 3) + distZ(:, :, :)*VoxelSize(3);

%% Smoothing of voxels
% Any remaining NaNs had to be inside, so we interpolate them (setting them to zeros creates artifacts!)
for iDim = 1:3
	IM(:, :, :, iDim) = xASL_im_ExtrapolateSmoothOverNaNs(IM(:, :, :, iDim), 0);
end

end
