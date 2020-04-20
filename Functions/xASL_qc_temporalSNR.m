function tSNR = xASL_qc_temporalSNR(pathIm4D,pathImTissueProb)
% Computes temporal SNR QC parameters in different compartments based on the timeseries data and tissue-probability masks
%
% FORMAT: tSNR = xASL_qc_temporalSNR(pathIm4D,pathImTissueProb)
%
% INPUT:
%   pathIm4D         - path to the 4D timeseries volume (REQUIRED)
%   pathImTissueProb - path to the 3D volume containing the tissue probability map (REQUIRED)
%                    - GM and WM paths needs to be given in a cell array of two strings, or a 2xN string
%
% OUTPUT:
%   tSNR     - struct containing the temporal SNR parameters for different compartments
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function computes several temporal SNR QC parameters as proposed in SPM Univariate Plus:
%            tSNR.tSNR_GM_Ratio      : mean GM signal / std GM over time 
%            tSNR.tSNR.tSNR_WM_Ratio : mean WM signal / std WM over time
%            tSNR.tSNR.tSNR_CSF_Ratio: mean CSF signal / std CSF over time
%            tSNR.tSNR_WMref_Ratio   : mean signal/std over time in eroded deep WM
%            tSNR.tSNR_GMWM_Ratio    : mean (GM+WM) signal / sqrt(std(GM+WM)^2+std(WMref)^2)
%            tSNR.tSNR_GMWM_WMref_Ratio: mean (GM+WM) signal / std WMref over time
%            tSNR.tSNR_Physio2Thermal_Ratio: sqrt((tSNR(GM+WM)/tSNR_GMWM_WMref_Ratio))^2-1)
%            tSNR.tSNR_Slope_Corr:
%            Differences to the SPM U+ suggestion: 
%            - eroded WM is used for estimating background noise
%            - Brainmask is determined in the same way as the structural anatQC,
%            - CSF is determined from the pGM&pWM maps;
%
% REFERENCES:
%              1) Thomas Liu (2016). Noise contributions to the fMRI signal: An overview NeuroImage, 343, 141-151
%                 http://dx.doi.org/10.1016/j.neuroimage.2016.09.008
%              2) Cesar Caballero-Gaudes and Richard C. Reynolds (2016). Methods For Cleaning The BOLD fMRI Signal. NeuroImage, 154,128-149
%              3) Lawrence Wald and Jonathan R Polimeni (2016). Impacting the effect of fMRI noise through 
%                 hardware and acquisition choices ??? Implications for controlling false positive rates. NeuroImage, 154,15-22
%              4) SPM Utility + toolbox. Cyril Pernet. https://osf.io/wn3h8/
%
% EXAMPLE: tSNR = xASL_qc_temporalSNR(x.P.Path_PWI4D,{x.P.Path_PVgm x.P.Path_PVwm});
% __________________________________
% Copyright (C) 2015-2019 ExploreASL

%% Input parameter admin

% Check if the input parameters are fine
if nargin ~= 2
	error('Need exactly two inputs.');
end

if isempty(pathIm4D)
	error('First input volume not given.');
end

if isempty(pathImTissueProb)
	error('Second input is empty.');
end

% Checks if two parameters are provided. And convert to cell if in array
if iscell(pathImTissueProb)
	if length(pathImTissueProb)~=2
		error(['Two tissue maps files expected, ' num2str(length(pathImTissueProb)) ' detected - check input file']);
	end
else
	if size(pathImTissueProb,1)~=2
		error(['Two tissue maps files expected, ' size(pathImTissueProb,1) ' detected - check input file']);
	end
	pathImTissueProb = mat2cell(pathImTissueProb,[1 1],size(pathImTissueProb,2));
end

% Make sure the Niftis are unpacked
xASL_adm_UnzipNifti(pathIm4D);
for iL=1:length(pathImTissueProb)
    xASL_adm_UnzipNifti(pathImTissueProb{iL});
end

%% Reading the data
% Read the time-series
% Check the names that they don't have ,1 at the end
vol4D = spm_vol(pathIm4D);
if  size(pathIm4D,1) == 1 && strcmp(pathIm4D(length(pathIm4D)-1:end),',1')
    pathIm4D = pathIm4D(1:length(pathIm4D)-2); % in case picked 4D put left ,1
end

nii4D     = xASL_io_ReadNifti(pathIm4D);

% Only run the QC when there are enough images
if iscell(vol4D)
	vol4D = cell2mat(vol4D); 
end

if size(vol4D,1) < 5
	warning('There are less than 5 images in your time series, skipping temporal QC');
	return;
end

% Load the PV maps
for m=1:2
	volTP(m) = spm_vol(pathImTissueProb{m});
end

if iscell(volTP)
	volTP = cell2mat(volTP);
end

% Check the dimensions
if any(volTP(1).dim~= volTP(2).dim)
    error('Tissue maps do not have the same dimentions');
end

if any(vol4D(1).dim~= volTP(1).dim)
	error('Dimension mismatch between the time series and tissue maps');
end

%% Compute tissue maps & masks
pGM         = spm_read_vols(volTP(1));
pWM         = spm_read_vols(volTP(2));
VoxelSize   = nii4D.hdr.pixdim(2:4);

% this filter needs to be twice as wide as we would with the anatomical high-resolution images (i.e. 8/voxelsize)
% because of the large voxel size of the functional images
BrainMask   = xASL_im_ndnanfilter(pGM+pWM,'gauss',double([16 16 16]./VoxelSize),0)>0.1;

pCSF        = max(0,BrainMask-pGM-pWM);

% Create mutually exclusive masks
GMmask = (pGM.*BrainMask)>0.5;
WMmask = (pWM.*BrainMask)>0.5;
CSFmask = (pCSF.*BrainMask)>0.5;

% Deep WM
WMeroded = xASL_im_DilateErodeFull(WMmask ,'erode',xASL_im_DilateErodeSphere(1));

%% Get temporal SNR (tSNR) within masks
tSNR = struct;

% GM SNR
DataN = xASL_Index2Data(vol4D, GMmask);
stdGM = xASL_stat_MeanNan(xASL_stat_StdNan(DataN,1,1));
tSNR.tSNR_GM_Ratio = xASL_stat_MeanNan(xASL_stat_MeanNan(DataN,1)) / stdGM;

% WM SNR
DataN = xASL_Index2Data(vol4D, WMmask);
tSNR.tSNR_WM_Ratio = xASL_stat_MeanNan(xASL_stat_MeanNan(DataN,1)) / xASL_stat_MeanNan(xASL_stat_StdNan(DataN,1,1));

% CSF SNR
DataN = xASL_Index2Data(vol4D, CSFmask);
tSNR.tSNR_CSF_Ratio = xASL_stat_MeanNan(xASL_stat_MeanNan(DataN,1)) / xASL_stat_MeanNan(xASL_stat_StdNan(DataN,1,1));

%% Get background tSNR (use WM as pseudo-background)
DataN = xASL_Index2Data(vol4D, WMeroded);
% tSNR.tSNR_WMref = xASL_stat_MeanNan(xASL_stat_MeanNan(DataN,1) ./ xASL_stat_StdNan(DataN,1,1) );
stdWMref = xASL_stat_MeanNan(xASL_stat_StdNan(DataN,1,1));
tSNR.tSNR_WMref_Ratio = xASL_stat_MeanNan(xASL_stat_MeanNan(DataN,1))/stdWMref;

%% tSNR of average brain tissue
DataN = xASL_Index2Data(vol4D, GMmask+WMmask);

meanGMWM = xASL_stat_MeanNan(DataN,1);
stdGMWM  = xASL_stat_MeanNan(xASL_stat_StdNan(DataN,1,1));
meanGMWM = xASL_stat_MeanNan(meanGMWM);
tSNR.tSNR_GMWM_Ratio = meanGMWM/sqrt(stdGMWM^2+stdWMref^2);
tSNR.tSNR_GMWM_WMref_Ratio = meanGMWM/stdWMref;

%% Physiothermal
tSNR.tSNR_Physio2Thermal_Ratio = abs(sqrt(((tSNR.tSNR_GMWM_Ratio/tSNR.tSNR_GMWM_WMref_Ratio)^2)-1));

%% Calculate SNR slope
IndexI = 1;
fprintf('Calculating SNR slope over voxels included for increasing brain probabilities:   ');
MaxP = 0.95;
MinP = 0.1;
StepP = 0.05;
nSteps = (MaxP-MinP)/StepP;
RangeP = xASL_round(MaxP:-StepP:MinP,2);
for p=RangeP
    CurrStep = find(RangeP==p);
    xASL_TrackProgress(CurrStep,nSteps);

    DataN = xASL_Index2Data(vol4D, pGM>=p);
    stdGM = xASL_stat_MeanNan(xASL_stat_StdNan(DataN,1,1));

    DataN = xASL_Index2Data(vol4D, pWM>p+pCSF>p);
    stdWMCSF = xASL_stat_MeanNan(xASL_stat_StdNan(DataN,1,1));

    ROI = pGM>p + pWM>p + pCSF>p;
    DataN = xASL_Index2Data(vol4D, ROI+(BrainMask~=1));

    ROIvalue(IndexI)= xASL_stat_MeanNan(DataN(:)) / sqrt(stdGMWM^2+stdWMref^2);
    ROIsize(IndexI) = sum(ROI(:));
    IndexI = IndexI + 1;
end

B = pinv([sqrt(ROIsize)' ones(18,1)])*ROIvalue';
tSNR.tSNR_Slope_Corr = B(1);
fprintf('\n');
end

%xASL_Index2Data Obtain data from coordinates as indicated by the imMask
function [DataN] = xASL_Index2Data(volIn,imMask)
    SizeIM = size(imMask);
    [x,y,z] = ind2sub(SizeIM,find(imMask));
    DataN = spm_get_data(volIn,[x y z]');
end
