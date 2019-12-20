function [DVARS] = xQC_DVARS(Bold)
%  xQC_DVARS: QC function to extract DVARS values
%
% FORMAT: [DVARS] = xQC_DVARS(NiiPath)
%
% INPUT:
%   Bold           - Can either be a 4D image or a path to the Bold.nii 
%   
% OUTPUT:
%   DVARS           - Vector of length T-1 containing DVARS values
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This Function computes DVARS QC values from a 4D image. 
%              D stans for Temporal Derivates of Timecourses
%              VARS is the RMS variance over Voxels 
%              DVARS indexes the rate of change across the entire brain at
%              each frame of data 
%
% 
% EXAMPLE: [DVARS] = xQC_DVARS(Path/To/Bold.nii)



BoldIm = xASL_io_Nifti2Im(Bold);

% Reshape into Timecourses Format
I = [size(BoldIm, 1),size(BoldIm, 2),size(BoldIm, 3)];
T = size(BoldIm, 4);
BoldTC = reshape(BoldIm,[prod(I), T]);

% Remove zeros and Nans
BoldTC = double(BoldTC);
zeros_idx  = find(sum(BoldTC,2)==0);
BoldTC(zeros_idx, :) = [];
nan_idx    = find(isnan(sum(BoldTC,2)));
BoldTC(nan_idx, :) = [];

% Intensity Normalization %
%No intensity normalization for the moment ... discuss
%MeanIm  = xASL_stat_MeanNan(BoldTC,2);
%Md = xASL_stat_MedianNan(MeanIm);

% Centre the data
MeanIm     =    xASL_stat_MeanNan(BoldTC,2);
dmeaner =    repmat(MeanIm,[1,T]);
BoldTC      =    BoldTC-dmeaner;

% DVARS
nVox = size(BoldTC, 1);
Deriv    = diff(BoldTC,1,2);
Deriv_2  = Deriv.^2;
MeanDer = xASL_stat_MeanNan(Deriv_2,1);
DVARS = sqrt(MeanDer);

end 







