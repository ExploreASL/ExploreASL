function [IM] = xASL_im_FillNaNs(InputPath, UseMethod, bQuality, VoxelSize, bZeroes)
%xASL_im_FillNaNs Fill NaNs in image
% FORMAT: xASL_im_FillNaNs(InputPath[, UseMethod, bQuality, VoxelSize, bZeroes])
%
% INPUT:
%   InputPath   - path to image, or image matrix (REQUIRED)
%                 NB: if this is a path, it is overwritten
%   UseMethod   - one of the following methods to fill the NaNs:
%                 1: extrapolate by smoothing (useful for any image, e.g.
%                    extrapolating M0 smooth biasfield) (DEFAULT)
%                 2: set to 0 (useful for e.g. DARTEL displacement fields,
%                    where 0 means no displacement)
%                 3: linearly extrapolate 
%                 4: fill with flow field identity map (specific to MNI 1.5mm space used by CAT12/DARTEL)
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
%   IM          - image matrix in which NaNs were filled
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

if UseMethod==1 % only method that needs VoxelSize
    if nargin<4 || isempty(VoxelSize)
        if ~isnumeric(InputPath) && xASL_exist(InputPath, 'file')
            nii = xASL_io_ReadNifti(InputPath);
            VoxelSize = nii.hdr.pixdim(2:4);
        elseif bQuality
                warning('Voxel size not specified, using default [1.5 1.5 1.5]');
                % only when high quality will we use this parameter,
                % so we dont give this warning for low quality
                VoxelSize = [1.5 1.5 1.5];
        else
                VoxelSize = [1.5 1.5 1.5];
        end
    end
end

if nargin<5 || isempty(bZeroes)
    if UseMethod==3
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


%% ------------------------------------------------------------------------------
%% 2. Replace NaNs using any of the methods
fprintf('%s', 'Cleaning up NaNs in image:  0%');
TotalDim = size(IM,4)*size(IM,5)*size(IM,6)*size(IM,7);
CurrentIt = 1;
for i4=1:size(IM,4)
    for i5=1:size(IM,5)
        for i6=1:size(IM,6)
            for i7=1:size(IM,7)
                xASL_TrackProgress(CurrentIt, TotalDim);

                switch UseMethod
                    case 1
                        % Use this for cases where NaNs need to get border values (e.g. extrapolating a masked image)
                        IM(:,:,:,i4,i5,i6,i7) = xASL_im_ExtrapolateSmoothOverNaNs(IM(:,:,:,i4,i5,i6,i7), bQuality, VoxelSize);
                    case 2
                        IM(isnan(IM)) = 0;
                        % This loops unnecessary but goes fast enough, keep in this loop for readability
                    case 3
                        % Use this for cases where NaNs need to expand linearly (e.g., extrapolating flowfields)
                        IM(:,:,:,i4,i5,i6,i7) = xASL_im_ExtrapolateLinearlyOverNaNs(IM(:,:,:,i4,i5,i6,i7));
                    otherwise
                        error('Invalid option');
                end

                CurrentIt = CurrentIt+1;
            end
        end
    end
end
fprintf('\n');


%% ------------------------------------------------------------------------------
%% 3. Save image
% Here we only output an image if the input path was a path and not an
% image matrix
if ~isnumeric(InputPath) && xASL_exist(InputPath, 'file')
    xASL_io_SaveNifti(InputPath, InputPath, IM, [], 0);
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
function [IM] = xASL_im_ExtrapolateLinearlyOverNaNs(IM, nEdgeVoxels)
%xASL_im_ExtrapolateLinearlyOverNaNs
% The input of this function is 3D (iterations of other dimensions are done
% in the main function above)
%
% IM          - input image
% nEdgeVoxels - number of voxels at the edges of the non-nan values that we check for weird interpolation values, 
% to be set to NaN and ignored in the extrapolation

%% 0. Admin
% Also catch Infs or imaginary numbers or whatever. Zeroes are already catched above in the main function.
IM(~isfinite(IM)) = NaN;

if nargin<2 || isempty(nEdgeVoxels)
    nEdgeVoxels = 0;
    % We can try this later, but it seems to work less well than doing the same based on the distance transform
end


%% 0.5 Fix interpolation edges
% The first layer(s) next to the NaNs are sometimes interpolated as well
% So here we try to estimate which voxels those are, by looking at the distance from the NaNs
% together with the transformation difference from a smoothed transformation field.

% First we get the distance transform
distMap = xASL_im_DistanceTransform(isnan(IM));

% Then we create a copy of the image that is smoothed (outside of NaNs)
IMsmooth = xASL_im_ndnanfilter(IM, 'gauss', [8 8 8], 1);
% We obtain an absolute difference image between the original transformation field, and the smoothed transformation field
% This difference will show us how "extreme" the transformations are.
deltaIm = abs(IM-IMsmooth);
% This absolute difference image, is now divided by the distance.
% Because we want to detect extreme transformation values that are close to the border (low distance)
% So this would result in a high relativeDelta
relativeDelta = deltaIm./distMap;

% The relativeDelta is Rician-ish distributed, so we use the median and MAD to obtain a threshold for the relativeDelta
medianRelativeDelta = xASL_stat_MedianNan(relativeDelta(:));
MADRelativeDelta = xASL_stat_MadNan(relativeDelta(:));
thresholdIs = medianRelativeDelta+3.5*MADRelativeDelta;

% We want to have a hard cutoff of distance from NaNs that we won't change
% This is 10% of the transformation field size (e.g, for 121*145*121 transformation field, this is 13 voxels)
sizeIm = size(IM);
if numel(sizeIm)<3
    warning('Something wrong with this flowfield');
end
averageSize = mean(sizeIm(1:3));
edgeDistance = round(0.1*averageSize);

% Set the voxels at the border where extreme values were detected to NaN as well
IM(relativeDelta>thresholdIs & distMap<edgeDistance) = NaN;

nX = size(IM,1);
nY = size(IM,2);
nZ = size(IM,3);

% Now we start extrapolating all NaNs — including those that we just created ourselves at the border — with robust non-NaN values at the border

for Iteration = 1:3 % Applies thrice along each X Y Z dimension to make sure all corners are filled
    if sum(isnan(IM(:)))>0
        %% 1. Reshape 3D to 2D
        % Along the three dimensions (or two dimensions in case a 2D image is on the input)
	    % Check each line of the image along the dimension and checks NaNs at the start and end of each line
	    % And interpolates them using the values in the middle
        for DimensionIs = 1:min(3,ndims(IM))
            % Takes a 3D image and concatenates to a 2D so that we can fill in the NaNs in each column
            switch(DimensionIs)
                case 1
                    imCol = reshape(IM, [nX,nY*nZ]);
                case 2
                    imCol = reshape(shiftdim(IM,1), [nY,nZ*nX]);
                case 3
                    imCol = reshape(IM, [nX*nY,nZ])';
            end
    
		    % Find columns with a maximum amount of NaNs.
            % Other NaNs will be dealt with by the xASL_im_ExtrapolateSmoothOverNaNs function below
            % 50% at first iteration
            % 10% at later iterations
		    if Iteration==1
			    indCol = sum(isnan(imCol),1)>0 & sum(~isnan(imCol),1)>(0.5*size(imCol,1));
		    else
			    indCol = sum(isnan(imCol),1)>0 & sum(~isnan(imCol),1)>(0.1*size(imCol,1));
		    end			
		    
            nCol = size(imCol,1);
            % Find columns starting or ending with nans
            indColLow  = indCol & isnan(imCol(1,:));
		    indColHigh = indCol & isnan(imCol(end,:));
    
            %% 2. Linearly extrapolate NaNs for columns starting with NaNs
            % For each column starting with Nans checks how the values are changing in the non-NaN voxels and linearly
		    % interpolate in the NaNs
            for ColumnIs = find(indColLow)
                % Take this column only
                imThisCol = imCol(:,ColumnIs);
                
                [diffIs, minInd, maxInd, thresholdIs] = xASL_im_Extrapolate_FixEdges(imThisCol, 0);

                % Now we check if we should remove edge voxels
                % We try this nEdgeVoxels times (default=5)
                nLayersProcessed = 0;
                while nLayersProcessed<nEdgeVoxels
                    if diffIs(end)>thresholdIs
                        imThisCol(minInd) = NaN; % remove the edge voxel
                        % If we have changed a first or last non-nan voxel, we rerun the difference determination 
                        [diffIs, minInd, maxInd, thresholdIs] = xASL_im_Extrapolate_FixEdges(imThisCol, 1);
                    end
                    nLayersProcessed = nLayersProcessed+1;
                end

                meanDiff = xASL_stat_MeanNan(diffIs);
                
                % Create a robust first value, from 6 voxels
                % This avoids extrapolating an extreme value at the edge
                voxels6 = imThisCol(minInd:minInd+5); % select the first 6 voxels
                voxels6 = voxels6(isfinite(voxels6)); % select those without NaNs
                robustMedian = median(voxels6);
                robustDiff = median(voxels6(1:end-1) - voxels6(2:end)); % mean difference
                numVoxels = length(voxels6)/2;
                robustValue = robustMedian + numVoxels*robustDiff;
                % Now take the mean of the actual first value and the robust one
                robustValue = mean([robustValue imThisCol(minInd)]);

                % Apply this
                imThisCol(1:minInd-1) = robustValue + (minInd-1:-1:1) * meanDiff;
                % Put this column back in the image
                imCol(:, ColumnIs) = imThisCol;
            end
            
            %% 3. Linearly extrapolate NaNs for columns ending with NaNs
            % For each column ending with NaNs check how the values are changing in the non-NaN voxels and linearly
		    % interpolate in the NaNs
            for ColumnIs = find(indColHigh)
                % Take this column only
                imThisCol = imCol(:,ColumnIs);
    
                % In this subfunction, we determine the difference within the current column
                % Which we can use for extrapolation, and to check if we need to exclude edge voxels
                % of the current column
                [diffIs, minInd, maxInd, thresholdIs] = xASL_im_Extrapolate_FixEdges(imThisCol, 1);

                % Now we check if we should remove edge voxels
                % We try this nEdgeVoxels times (default=5)
                nLayersProcessed = 0;
                while nLayersProcessed<nEdgeVoxels
                    % 
                    if diffIs(end)>thresholdIs
                        imThisCol(maxInd) = NaN; % remove the edge voxel
                        % If we have changed a first or last non-nan voxel, we rerun the difference determination 
                        [diffIs, minInd, maxInd, thresholdIs] = xASL_im_Extrapolate_FixEdges(imThisCol, 1);
                    end
                    nLayersProcessed = nLayersProcessed+1;
                end

                meanDiff = xASL_stat_MeanNan(diffIs);

                % Create a robust last value, from 6 voxels
                % This avoids extrapolating an extreme value at the edge
                voxels6 = imThisCol(maxInd-5:maxInd); % select the last 6 voxels
                voxels6 = voxels6(isfinite(voxels6)); % select those without NaNs
                robustMedian = median(voxels6);
                robustDiff = median(voxels6(2:end) - voxels6(1:end-1)); % mean difference
                numVoxels = length(voxels6)/2;
                robustValue = robustMedian + numVoxels*robustDiff;
                % Now take the mean of the actual last value and the robust one
                robustValue = mean([robustValue imThisCol(maxInd)]);

                % Apply this
                imThisCol(maxInd+1:end) = robustValue + (1:nCol-maxInd) * meanDiff;
                % Put this column back in the image
                imCol(:, ColumnIs) = imThisCol;
            end
    
            %% 4. Reshape 2D back to 3D
            switch(DimensionIs)
                case 1
                    IM = reshape(imCol, [nX,nY,nZ]);
                case 2
				    % For 2D matrix and 3D matrix, the reshaping back is done differently
				    if ndims(IM)<3
					    IM = shiftdim(reshape(imCol,[nY,nX]),1);
				    else
					    IM = shiftdim(reshape(imCol,[nY,nZ,nX]),2);
				    end
                case 3
                    IM = reshape(imCol', [nX,nY,nZ]);
            end % switch
        end % for DimensionIs
    end % sum(isnan(IM))
end

% Any remaining NaNs had to be inside, so we interpolate them (setting them
% to zeros creates artifacts!)
IM = xASL_im_ExtrapolateSmoothOverNaNs(IM, 0);

end


function [diffIs, minInd, maxInd, thresholdIs] = xASL_im_Extrapolate_FixEdges(imThisCol, bSide)
%xASL_im_Extrapolate_FixEdges Check if the edges of the non-NaN voxels within a column are outliers
% from the distribution of differences between voxels. If they are outliers, they are probably interpolation 
% errors from the SPM software, and we want to ignore them.
%
% input
% imThisCol - vector of current column (REQUIRED)
% bSide     - boolean specifying if we check the start/first non-NaN voxel (0)
% %           or the end/last non-NaN voxel (1). (REQUIRED)


minInd = find(~isnan(imThisCol), 1, 'first'); % Get the index of the first non-value
maxInd = find(~isnan(imThisCol), 1, 'last'); % Get the index of the last non-value
midInd = floor(minInd+(maxInd-minInd)/2); % Get the index of the middle voxel

% Mean difference, from the middle index to the index of the last value
% Calculated by their difference with their neighboring voxel            

if ~bSide % start/first
    diffIs = imThisCol(minInd:midInd) - imThisCol(minInd+1:midInd+1);
else
    diffIs = imThisCol(midInd:maxInd) - imThisCol(midInd-1:maxInd-1);
end

% Remove outlier voxels near the edges, due to interpolation
% When a voxel near the edge has 1.96*SD larger diff than the mean of absolute diff, we set it to NaN
% For a maximum of maxEdgeVoxels

absDiff = abs(diffIs);
meanAbsDiff = xASL_stat_MeanNan(absDiff);
stdDiff = xASL_stat_StdNan(diffIs);
thresholdIs = meanAbsDiff+3.5*stdDiff;

end