function [NotOutliers, iOutliers, RMS] = xASL_stat_RobustMean(IM, ParameterFunction, imMask)
% Submodule of ExploreASL Structural module, that obtains volumes from the tissue segmentations
% (& FLAIR WMH segmentations if they exist)
%
% FORMAT: [NoOutliers, iOutliers, RMS] = xASL_stat_RobustMean(IM, ParameterFunction)
%
% INPUT:
%   IM                 - input images, should be masked with xASL_im_IM2Column: vector image (REQUIRED)
%   ParameterFunction  - parametric function to use for defining deviation
%                        of an image: options:
%                                     SoS - sum of squared errors (DEFAULT)
%                                     AI  - average relative asymmetry index
%   imMask             - for parametric maps (e.g. CBF, ATT, Tex, DCE Ktrans) it can be more useful to only inspect within the GM and WM
%                        as poor fits in the CSF or outside the brain are less interesting (OPTIONAL, DEFAULT=whole image)
%
% OUTPUT:
%   NotOutliers         - vector, true for images that were not outliers
%   iOutliers           - indices of images that were outliers
%   RMS                 - vector of numerical values, RMS of difference of individual image with group-average
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function detects outlier images, that can be used to create
%              a robust average, e.g. for template or biasfield creation. This is based either on the sum-of-squares
%              with the mean image (SoS), or on the average relative asymmetry index (AI). Images that are 
%              median+/-3 mad off are defined as outliers. MAD = median/mean absolute difference
%
% EXAMPLE: NotOutliers = xASL_stat_RobustMean(IM);
% __________________________________
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________



%% Admin
% if size(IM,2)<16 % only do the outlier detection with sufficiently large datasets
%     NotOutliers = ones(size(IM,2),1);
%     iOutliers = []; % empty
%     RMS = NaN;
%     fprintf('Outlier exclusion skipped, too small dataset\n');
%     return;
% end
if nargin<3 || isempty(imMask)
    imMask = true([size(IM, 1), 1]);
end

if nargin<2 || isempty(ParameterFunction)
    ParameterFunction = 'SoS';
elseif isempty(regexpi(ParameterFunction, '^(SoS|AI)$'))
    warning(['Unknown ParameterFunction: ' ParameterFunction ', using SoS']);
    ParameterFunction = 'SoS';
end
if size(IM,2)>size(IM,1) || ndims(IM)>2
    warning('Input IM has incorrect dimensions');
end

nScans = size(IM, 2);

%% Mask the image
if size(imMask, 1)~= size(IM, 1)
    error('Image has a different size than the mask, different spaces?');
elseif ndims(imMask)>2 || size(imMask, 2)~=1
    error('Mask should have a single dimension');
end

IM = IM(imMask, :);


%% Compute median, MAD, & deviations


% Create template image, to compare with
fprintf('%s\n',['QC: detecting outliers for n=' num2str(nScans)]);

MedianIM = repmat(xASL_stat_MedianNan(IM, 2), [1 nScans]);

if strcmpi(ParameterFunction,'SoS')
     DiffIm = (IM - MedianIM).^2;
     Deviation = xASL_stat_MeanNan(DiffIm, 1); % gives deviation sum per image, higher is worse quality
     Deviation = sqrt(Deviation);
     RMS = Deviation;
elseif strcmpi(ParameterFunction,'AI')
     DiffIm = abs(IM - repmat(IMtemp,[1 nScans])) ./ (0.5.*(IM + repmat(IMtemp,[1 nScans]))); % weighted SoS, AI
     Deviation = xASL_stat_MeanNan(DiffIm, 1); % gives deviation sum per image, higher is worse quality
end

NaNmask = isfinite(Deviation);
MedianDeviation = median(Deviation(NaNmask)); % average deviation from average image
MadDeviation = median(abs(Deviation(NaNmask) - MedianDeviation)); % Mean Absolute Difference (MAD)


%% Compute threshold & provide indices for those that are not outliers (i.e. not above threshold)
ThresholdDeviation = MedianDeviation+3.*MadDeviation;
NotOutliers = ~(Deviation>ThresholdDeviation)';
iOutliers = find(Deviation>ThresholdDeviation)';

if iOutliers>0
    fprintf(['Detected ' num2str(numel(iOutliers)) ' outliers\n']);
end


end