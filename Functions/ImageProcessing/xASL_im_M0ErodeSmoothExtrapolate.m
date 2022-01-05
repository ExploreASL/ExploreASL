function [ImOut, VisualQC] = xASL_im_M0ErodeSmoothExtrapolate(ImIn, DirOutput, NameOutput, pvGM, pvWM, brainCentralityMap, LowThreshold)
%xASL_im_M0ErodeSmoothExtrapolate M0 image processing
%
% FORMAT: [ImOut, VisualQC] = xASL_im_M0ErodeSmoothExtrapolate(ImIn, DirOutput, NameOutput, pvGM, pvWM, brainCentralityMap[, LowThreshold])
%
% INPUT:
%   ImIn - unprocessed M0 image, as 3D image or path (REQUIRED)
%   DirOutput - sring path to output folder, used to be x.D.M0regASLdir (REQUIRED)
%   NameOutput - string for filename in ['M0_im_proc_' NameOutput '.jpg'] (REQUIRED)
%   pvGM  - unprocessed pvGM image (3D image or path, same space as M0 image (REQUIRED)
%   pvWM  - unprocessed pvWM image (3D image or path, same space as M0 image (REQUIRED)
%   brainCentralityMap - brain mask probability map multiplied by centrality map, as 3D image or path (REQUIRED)
%   LowThreshold - numerical value between 0 and 1 for percentile sorted images values that will define the inclusion mask (mask>LowThreshold)
%                   (OPTIONAL, DEFAULT=0.7)
%
% OUTPUT:
%   ImOut    - processed M0 image
%   VisualQC - visualization of the M0 processing steps for quality control
% OUTPUT FILES: Visual QC image of M0 processing
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function erodes, smooths & extrapolates M0 in MNI standard space (tested at 1.5x1.5x1.5mm as used in ExploreASL).
%               It assumes that the M0 image is in standard space & that the GM & WM probability maps
%               are aligned. Here, we mask the M0, to remove high CSF signal and low extracranial signal,
%               enabling us to smooth the image without inserting wrong signal. See also
%               the ExploreASL manuscript for a more extensive explanation. 
%
%               A visual QC figure is created, showing the M0 image processing steps for a single transversal slice (slice 53 in `1.5 mm` MNI standard space)
%               `OutputFile = fullfile(DirOutput,['M0_im_proc_' NameOutput '.jpg']);`
%               The original M0 image (a) is masked with a (pGM+pWM)>50% mask (b)
%               eroded with a two-voxel sphere to limit the influence of the ventricular and extracranial signal (c)
%               and thresholded to exclude significantly high (i.e. median + 3*mean absolute deviation (MAD)) border region values (d)
%               This masked M0 image is smoothed with a 16x16x16 mm full- width-half-maximum Gaussian filter (Mutsaerts et al., 2018) (e)
%               after which the signal at the border is smoothly extrapolated until the full image is filled (f).
%               Whereas the masking avoids mixing with cerebrospinal fluid or extracranial signal, the extrapolation avoids M0 division artifacts
%
%               This function performs the following steps:
%
%               1.  Mask 1: structural mask
%               2.  Mask 2: intensity-based mask to remove extracranial signal
%               2.a Multiplication with brain mask centrality map
%               2.b Mask the "dummy image" from 2.a based on intensities
%               2.c Enforce to include the eroded brain in the mask
%               3.  Mask 3: Erode the combined masks
%               4.  Mask 4: Determine any odd borders
%               5.  Smoothing
%               6.  Extrapolating only
%               7.  Scale back to the GM M0
%               8.  Print visual QC figure
%
% EXAMPLE: as used in ExploreASL (xASL_wrp_ProcessM0): [ImOut] = xASL_im_M0ErodeSmoothExtrapolate(x.P.Pop_Path_M0, x.D.M0regASLdir, x.P.SubjectID, x.P.Pop_Path_rc1T1, x.P.Pop_Path_rc2T1, path_BrainCentralityMap)
% __________________________________
% Copyright (C) 2015-2021 ExploreASL


%% ------------------------------------------------------------------------------------------
%% 0. Admin
if nargin < 1 || isempty(ImIn)
    error('Please specify M0 image to be processed');
elseif nargin < 2 || isempty(DirOutput)
	error('Please specify output folder');
elseif nargin < 3 || isempty(NameOutput)
	error('Please specify the output name');
elseif nargin < 4 || isempty(pvGM)
	error('Please specify the pvGM');
elseif nargin < 5 || isempty(pvWM)
	error('Please specify the pvWM');
elseif nargin < 6 || isempty(brainCentralityMap)
	error('Please specify the brainCentralityMap');
elseif nargin < 7 || isempty(LowThreshold)
    LowThreshold = 0.7;
elseif ~isnumeric(LowThreshold)
    error('Invalid LowThreshold input parameter, should be numerical');
elseif LowThreshold<0 || LowThreshold>1
    error('Invalid LowThreshold value, should be between 0 and 1');
end

% Load the M0 image
ImIn = xASL_io_Nifti2Im(ImIn);


%% ------------------------------------------------------------------------------------------
%% 1. Mask 1: structural mask
% Here, we load the partial volume maps for the first mask
GMim = xASL_io_Nifti2Im(pvGM);
WMim = xASL_io_Nifti2Im(pvWM);

if xASL_stat_SumNan(GMim(:))==0
    error('Empty GM partial volume map, cannot process the M0');
elseif xASL_stat_SumNan(WMim(:))==0
    error('Empty WM partial volume map, cannot process the M0');
else
    GMmask = GMim>0.7;
    Mask1 = (GMim+WMim)>0.5;
    fprintf('\n%s', 'Processing M0 image: ');
end


%% ------------------------------------------------------------------------------------------
%% 2. Mask 2: intensity-based mask to remove extracranial signal
% 2.a Multiplication with brain mask centrality map
% We multiply the M0 image with the brain mask centrality map, to
% decrease the values outside the brain in intensity, so the
% intensity-based masking will be more biased towards the outside of the
% brain. This helps in cases with low intensities inside the brain e.g.,
% due to a biasfield. We call the resulting map here a "dummy image"
dummyImage = ImIn.*xASL_io_Nifti2Im(brainCentralityMap);

% Get NaNs
SortInt = sort(dummyImage(:));
SortInt = SortInt(~isnan(SortInt));
if isempty(SortInt)
    error('M0 image only contains NaN values, cannot process M0');
end

% 2.b Mask the "dummy image" from 2.a based on intensities
%   *   The high intensity masking removes residual high CSF values (for
%       T2-weighted M0)
%   *   The low intensity masking does an extra round of brainmasking,
%       which is useful with e.g., susceptibility artifacts,
%       in which cases the ASL brainmask should be smaller than the tissue
%       brainmask created in 1.
fprintf('%s', 'intensity masking, ');
% Lower threshold
ThresholdN = SortInt(round(LowThreshold*length(SortInt)));

% Remove peak signal as well
ThresholdN2 = SortInt(round(0.999*length(SortInt)));

% Combine these masks
Mask2 = (dummyImage>ThresholdN) & (dummyImage<ThresholdN2);
Mask2 = Mask1 & Mask2;

% 2.c Enforce to include the eroded brain in the mask
% Sometimes, the low intensity masking in 2.b can still remove brain
% regions. Here, we force that a strong eroded tissue brainmask is still
% restored (we don't expect susceptibility artifacts in the brainmask)
Mask1Eroded = xASL_im_DilateErodeFull(Mask1, 'erode', xASL_im_DilateErodeSphere(4));
% But it still has to fall within the FoV, so we mask it with NaNs
Mask1Eroded(isnan(ImIn)) = 0;
Mask2(Mask1Eroded) = 1;


%% ------------------------------------------------------------------------------------------
%% 3. Mask 3: Erode the combined masks
% We erode the resulting mask, as misalignment of ASL with the tissue
% masks could still result in suboptimal masking. To improve the masking here, we erode with a 2-voxel sphere. 
fprintf('%s', 'eroding, ');
Mask3 = xASL_im_DilateErodeFull(Mask2, 'erode', xASL_im_DilateErodeSphere(2));


%% ------------------------------------------------------------------------------------------
%% 4. Mask 4: Determine any odd borders
% We remove sharp edges, that still have local intensity outliers

% First get median inside the mask
ValuesM0mask = ImIn(Mask3 & isfinite(ImIn));
MedianN = median(ValuesM0mask);
% Fill voxels outside mask with this median
TempIM = ImIn;
TempIM(~Mask3) = MedianN;

% Smooth ImOut, to check the voxels that change most by smoothing
% (especially those with low values)
ImOutSmooth = xASL_im_ndnanfilter(TempIM, 'gauss', double([5 5 5]), 0);
DiffIM = abs(ImOutSmooth-TempIM);

DiffIMvalues = DiffIM(Mask3 & isfinite(DiffIM));

MedianN = median(DiffIMvalues);
MadN = xASL_stat_MadNan(DiffIMvalues);
HiThresh = MedianN+4.5*MadN;
% Here, we may remove too much, but this will be smoothed & extrapolated, anyway.
Mask4 = Mask3 & ~(DiffIM>HiThresh);


%% ------------------------------------------------------------------------------------------
%% 5. Smoothing
% Smooth M0 map before division
% Currently, same smoothing kernel for 2D & 3D, to remove e.g. Gibbs
% ringing artifact in 3D data. And smoothing sums quadratically

% This should be considerable large smoothing,
% because of the noise that would be introduced by
% voxel-wise division. See Beaumont's thesis, chapter 4

% Initial large smoothing for smooth brain image

fprintf('%s', 'smoothing, ');
ImOut = ImIn.*Mask4;
ImOut(ImOut==0) = NaN;

% smoothing with interpolation
% The bigger the kernel size, the smoother the biasfield, but the
% more artifacts will be filtered into the data

VoxelSize = [1.5 1.5 1.5];
ImOut = xASL_im_ndnanfilter(ImOut,'gauss',double([16 16 16]./VoxelSize),0);
Im5 = ImOut;
ImOut = xASL_im_ndnanfilter(ImOut,'gauss',double([16 16 16]./VoxelSize),0);


%% ------------------------------------------------------------------------------------------
%% 6. Extrapolating only
% We fill the residual NaNs (outside the mask) of the FoV
% to prevent ASL/M0 division artifacts
ImOut = xASL_im_FillNaNs(ImOut, 1, 1, VoxelSize);


%% ------------------------------------------------------------------------------------------
%% 7. Scale back to the GM M0
% Our relatively strict masking has mostly excluded the GM M0 and
% substituted this with WM M0 values as these are more robust against
% misalignment. To improve quantitative accuracy, we scale the M0 values
% back to the median of GM M0.
% Note that this is optimal for GM ASL but may be suboptimal for WM ASL
OldGMM0 = ImIn(GMmask & isfinite(ImIn));
OldGMM0 = median(OldGMM0);

NewGMM0 = ImOut(GMmask & isfinite(ImOut));
NewGMM0 = median(NewGMM0);

RatioN = OldGMM0/NewGMM0;
ImOut = ImOut.*RatioN;


%% ------------------------------------------------------------------------------------------
%% 8. Print visual QC figure
S2S = 53; % slice to show
VisualQC = [xASL_im_rotate(ImIn(:,:,S2S),90) xASL_im_rotate(ImIn(:,:,S2S).*Mask2(:,:,S2S),90) xASL_im_rotate(ImIn(:,:,S2S).*Mask3(:,:,S2S),90) ; xASL_im_rotate(ImIn(:,:,S2S).*Mask4(:,:,S2S),90) xASL_im_rotate(Im5(:,:,S2S),90) xASL_im_rotate(ImOut(:,:,S2S),90)];

xASL_adm_CreateDir(DirOutput);
OutputFile = fullfile(DirOutput,['M0_im_proc_' NameOutput '.jpg']);
fprintf('%s\n',['Please check visual QC: ' OutputFile]);
xASL_vis_Imwrite(VisualQC, OutputFile);


end