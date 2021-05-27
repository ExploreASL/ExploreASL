function [NotOutliers, iOutliers] = xASL_stat_RobustMean(IM, ParameterFunction)
% Submodule of ExploreASL Structural module, that obtains volumes from the tissue segmentations
% (& FLAIR WMH segmentations if they exist)
%
% FORMAT: [NoOutliers, iOutliers, ThresholdDeviation] = xASL_stat_RobustMean(IM, ParameterFunction)
%
% INPUT:
%   IM                 - input images, should be masked with xASL_im_IM2Column (REQUIRED)
%   ParameterFunction  - parametric function to use for defining deviation
%                        of an image: options:
%                                     SoS - sum of squared errors (DEFAULT)
%                                     AI  - average relative asymmetry index
%
% OUTPUT:
%   NotOutliers         - vector, true for images that were not outliers
%   iOutliers           - indices of images that were outliers
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function detects outlier images, that can be used to create
%              a robust average, e.g. for template or biasfield creation. This is based either on the sum-of-squares
%              with the mean image (SoS), or on the average relative asymmetry index (AI). Images that are 
%              median+/-3 mad off are defined as outliers. MAD = median/mean absolute difference
%
% EXAMPLE: NotOutliers = xASL_stat_RobustMean(IM);
% __________________________________
% Copyright 2015-2019 ExploreASL


%% Admin
if size(IM,2)<16 % only do the outlier detection with too small datasets
    NotOutliers = ones(size(IM,2),1);
    iOutliers = []; % empty
    ThresholdDeviation = Inf;
    RobustMean = xASL_stat_MeanNan(IM,2);
    fprintf('Outlier exclusion skipped, too small dataset\n');
    return;
end
if nargin<2 || isempty(ParameterFunction)
    ParameterFunction = 'SoS';
elseif isempty(regexp(ParameterFunction, '^(SoS|AI)$'))
    warning(['Unknown ParameterFunction: ' ParameterFunction ', using SoS']);
    ParameterFunction = 'SoS';
end
if size(IM,2)>size(IM,1)
    warning('Input IM has incorrect dimensions');
end


%% THIS IS PROBABLY REDUNDANT CODE
% if prod(size(IM))==1 % assume this is a memory mapping file
%     IM = shiftdim(IM.Data.data,1);
% end
% 
% IM2 = xASL_im_IM2Column(IM,x.utils.WBmask); % Make sure to convert to column if too large

%% Compute median, MAD, & deviations

Size4 = size(IM,2);

% Create template image, to compare with
fprintf('%s\n',['Detecting outliers for n=' num2str(Size4)]);

MedianIM = repmat(xASL_stat_MedianNan(IM,2),[1 Size4]);

if strcmp(ParameterFunction,'SoS')
     ValueMask = abs(IM - MedianIM);
elseif strcmp(ParameterFunction,'AI')
     ValueMask = abs(IM - repmat(IMtemp,[1 Size4])) ./ (0.5.*(IM + repmat(IMtemp,[1 Size4]))); % weighted SoS, AI
end        

Deviation = xASL_stat_MeanNan(ValueMask,1); % gives deviation sum per image, higher is worse quality
% so the higher deviation, the more outlier the image is
NaNmask = isfinite(Deviation);
MedianDeviation = median(Deviation(NaNmask)); % average deviation from average image
MadDeviation = median(abs(Deviation(NaNmask) - MedianDeviation)); % Mean Absolute Difference (MAD)

%% Compute threshold & provide indices for those that are not outliers (i.e. not above threshold)
ThresholdDeviation = MedianDeviation+3.*MadDeviation;
NotOutliers = ~(Deviation>ThresholdDeviation)';
iOutliers = find(Deviation>ThresholdDeviation)';

fprintf(['Detected ' num2str(numel(iOutliers)) ' outliers\n']);


end