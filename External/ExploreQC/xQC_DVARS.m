function [DVARS] = xQC_DVARS(NiiPath)
% xQC_DVARS: QC function to extract DVARS values
%
% FORMAT: [DVARS] = xQC_DVARS(NiiPath)
%
% INPUT:
%   NiiPath         - Can either be a 4D image or a path to the Bold.nii 
%   
% OUTPUT:
%   DVARS           - Vector of length size(IM4D)-1 containing DVARS values
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This Function computes DVARS QC values from a 4D image (IM4D). 
%              D stands for Temporal Derivates of Timecourses
%              VARS is the RMS variance over Voxels 
%              DVARS indexes the rate of change across the entire brain at
%              each frame of data.
%              So in other words, the DVARS is the spatial RMS of
%              differenced timecourse, or simply the mean voxel-wise
%              motion.
%
% 
% EXAMPLE: DVARS = xQC_DVARS('Path/To/Bold.nii');



IM4D = xASL_io_Nifti2Im(NiiPath);

% Reshape image into column
NumElements = numel(IM4D(:,:,:,1,1,1,1,1));
DataColumn = reshape(IM4D,[NumElements, size(IM4D, 4)]);

% Remove zeros and Nans
IndicesZeros  = find(sum(DataColumn,2)==0);
DataColumn(IndicesZeros, :) = [];
NaNIndices    = find(isnan(sum(DataColumn,2)));
DataColumn(NaNIndices, :) = [];

% Intensity Normalization %
%No intensity normalization for the moment ... discuss
%MeanIm  = xASL_stat_MeanNan(BoldTC,2);
%Md = xASL_stat_MedianNan(MeanIm);

% Centre the data
MeanIm = xASL_stat_MeanNan(DataColumn,2);
DeMeaner = repmat(MeanIm,[1, size(IM4D, 4)]);
DataColumn = DataColumn-DeMeaner;

% DVARS
DerivativeIs = DataColumn(:,2:end)-DataColumn(:,1:end-1); % calculate difference
DerivativeIs = DerivativeIs.^2; % square
DerivativeIs = xASL_stat_MeanNan(DerivativeIs,1); % spatial mean
DerivativeIs = sqrt(DerivativeIs); % root
DVARS = squeeze(DerivativeIs); % manage dimensionality

end 







