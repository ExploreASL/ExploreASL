function [ImOut, VisualQC] = xASL_im_M0ErodeSmoothExtrapolate(ImIn, DirOutput, NameOutput, pvGM, pvWM, Quality)

%xASL_im_M0ErodeSmoothExtrapolate M0 image processing
%
% FORMAT: [ImOut] = xASL_im_M0ErodeSmoothExtrapolate(ImIn, x)
%
% INPUT:
%   ImIn - unprocessed M0 image (3D image or path, REQUIRED)
%   DirOutput - sring path to output folder, used to be x.D.M0regASLdir (REQUIRED)
%   NameOutput - string for filename in ['M0_im_proc_' NameOutput '.jpg'], used to be x.P.SubjectID (REQUIRED)
%   pvGM  - unprocessed pvGM image (3D image or path, same space as M0 image, REQUIRED, used to be x.P.Pop_Path_rc1T1)
%   pvWM  - unprocessed pvWM image (3D image or path, same space as M0 image, REQUIRED, used to be x.P.Pop_Path_rc2T1)
%   Quality - boolean for high (1) or low (0) processing quality (used to be x.settings.Quality, OPTIONAL, DEFAULT = true)
%
% OUTPUT:
%   ImOut    - processed M0 image
%   VisualQC - visualization of the M0 processing steps for quality control
% OUTPUT FILES: Visual QC image of M0 processing
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function erodes, smooths & extrapolates M0 in standard space.
%               It assumes that the M0 image is in standard space & that the GM & WM probability maps
%               are aligned. Here, we mask the M0, to remove high CSF signal and low extracranial signal,
%               enabling us to smooth the image without inserting wrong signal. See also
%               the ExploreASL manuscript for a more extensive explanation. This function
%               performs the following steps:
%
%               1. Mask: Load segmentations, create structural mask
%               2. Mask: Create intensity-based mask to remove extracranial signal
%               3. Mask: Erode the combined masks
%               4. Mask: Determine any odd borders
%               5. Smoothing
%               6. Extrapolating only
%               7. Scale back to the GM M0
%               8. Print visual QC figure
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
% EXAMPLE: [ImOut] = xASL_im_M0ErodeSmoothExtrapolate(ImIn, x)
% __________________________________
% Copyright (C) 2015-2021 ExploreASL


%% ------------------------------------------------------------------------------------------
%% Admin)
if nargin<6 || isempty(Quality)
    fprintf('%s\n', 'Quality parameter missing, defaulting to normal quality');
    Quality = 1;
end

if nargin<3 || isempty(NameOutput)
    error('Please specify the output name');
elseif nargin<2 || isempty(DirOutput)
    error('Please specify output folder');
elseif nargin<1 || isempty(ImIn)
    error('Please specify M0 image to be processed');
end

% Initialize the defaults
ImIn = xASL_io_Nifti2Im(ImIn);

GMmask = zeros(size(ImIn));
Mask1 = zeros(size(ImIn));
Mask2 = zeros(size(ImIn));
ExistpGMpWM = 0;


%% ------------------------------------------------------------------------------------------
%% Mask 1) Load segmentations, create structural mask

if nargin<5 || isempty(pvWM) || isempty(pvGM)
    fprintf('%s\n','Masking M0 with intensity-based mask only, structural (pGM+pWM) files missing');
else
    GMim = xASL_io_Nifti2Im(pvGM);
    WMim = xASL_io_Nifti2Im(pvWM);
    
    if xASL_stat_SumNan(GMim(:))==0
        warning('Empty GM image');
    elseif xASL_stat_SumNan(WMim(:))==0
        warning('Empty WM image');
    else
        GMmask = GMim>0.7;
        Mask1 = (GMim+WMim)>0.5;
        ExistpGMpWM = 1;
    end
end


%% ------------------------------------------------------------------------------------------
%% Mask 2) Create intensity-based mask to remove extracranial signal
% Get NaNs
SortInt = sort(ImIn(:));
SortInt = SortInt(~isnan(SortInt));

% Create masks
if ~isempty(SortInt)
	% These thresholds are relatively conservative, otherwise this won't
	% work in cases of biasfields
    ThresholdN = SortInt(round(0.7*length(SortInt)));

	% Remove peak signal as well
	ThresholdN2 = SortInt(round(0.999*length(SortInt)));

	% Combine these masks
	if ExistpGMpWM
		Mask2 = Mask1 & (ImIn>ThresholdN) & (ImIn<ThresholdN2);
	else
		Mask2 = (ImIn>ThresholdN) & (ImIn<ThresholdN2);
		GMmask = ImIn>SortInt(round(0.85*length(SortInt))) & ImIn<SortInt(round(0.95*length(SortInt)));
	end
else
	warning('M0 image only contains NaN values...');
end

%% ------------------------------------------------------------------------------------------
%% Mask 3) Erode the combined masks
fprintf('%s\n','Erode M0 mask with 2-voxel sphere');
Mask3 = xASL_im_DilateErodeFull(Mask2, 'erode', xASL_im_DilateErodeSphere(2));

%% ------------------------------------------------------------------------------------------
%% Mask 4) Determine any odd borders
fprintf('%s\n','Identify & remove non-smooth values at the border of the M0 mask');

% First get median inside the mask
ValuesM0mask = ImIn(Mask3 & isfinite(ImIn));
MedianN = median(ValuesM0mask);
% Fill voxels outside mask with this median
TempIM = ImIn;
TempIM(~Mask3) = MedianN;

% Smooth ImOut, to check the voxels that change most by smoothing
% (especially those with low values)
ImOutSmooth = xASL_im_ndnanfilter(TempIM,'gauss',double([5 5 5]),0);
DiffIM = abs(ImOutSmooth-TempIM);

DiffIMvalues = DiffIM(Mask3 & isfinite(DiffIM));

MedianN = median(DiffIMvalues);
MadN = xASL_stat_MadNan(DiffIMvalues);
HiThresh = MedianN+4.5*MadN;
% here we remove too many parts, but this is not bad, as it will be smoothed & extrapolated, anyway.
Mask4 = Mask3 & ~(DiffIM>HiThresh);

%% ------------------------------------------------------------------------------------------
%% 5) Smoothing
% Smooth M0 map before division
% Currently, same smoothing kernel for 2D & 3D, to remove e.g. Gibbs
% ringing artifact in 3D data. And smoothing sums quadratically

% This should be considerable large smoothing,
% because of the noise that would be introduced by
% voxel-wise division. See Beaumont's thesis, chapter 4

% Initial large smoothing for smooth brain image

MaxIt = 24; % usually around this value (used for counting, not for limitation
                                % the number of iterations)

fprintf('Mask M0 with this mask & smooth...  ');
ImOut = ImIn.*Mask4;
ImOut(ImOut==0) = NaN;

% smoothing with interpolation
% The bigger the kernel size, the smoother the biasfield, but the
% more artifacts will be filtered into the data

VoxelSize = [1.5 1.5 1.5];

ImOut = xASL_im_ndnanfilter(ImOut,'gauss',double([16 16 16]./VoxelSize),0);
xASL_TrackProgress(1,MaxIt);
Im5 = ImOut;

if Quality
    ImOut = xASL_im_ndnanfilter(ImOut,'gauss',double([16 16 16]./VoxelSize),0);
else % in case of low quality, we leave this to the xASL_im_FillNaNs below,
     % where a smaller kernel will go much faster for low quality 
end
xASL_TrackProgress(2, MaxIt);

%% ------------------------------------------------------------------------------------------
%% 6) Extrapolating only
% Here we fill the residual NaNs (outside the mask) of the FoV
% to prevent ASL/M0 division artifacts
ImOut = xASL_im_FillNaNs(ImOut, 1, Quality, VoxelSize);


%% ------------------------------------------------------------------------------------------
%% 7) Scale back to the GM M0
fprintf('\n%s\n','Rescale the smooth biasfield GM M0 values back to non-smooth GM M0 values');
OldGMM0 = ImIn(GMmask & isfinite(ImIn));
OldGMM0 = median(OldGMM0);

NewGMM0 = ImOut(GMmask & isfinite(ImOut));
NewGMM0 = median(NewGMM0);

RatioN = OldGMM0/NewGMM0;
ImOut = ImOut.*RatioN;

%% ------------------------------------------------------------------------------------------
%% 8) Print visual QC figure

S2S = 53; % slice to show
VisualQC = [xASL_im_rotate(ImIn(:,:,S2S),90) xASL_im_rotate(ImIn(:,:,S2S).*Mask2(:,:,S2S),90) xASL_im_rotate(ImIn(:,:,S2S).*Mask3(:,:,S2S),90) ; xASL_im_rotate(ImIn(:,:,S2S).*Mask4(:,:,S2S),90) xASL_im_rotate(Im5(:,:,S2S),90) xASL_im_rotate(ImOut(:,:,S2S),90)];

xASL_adm_CreateDir(DirOutput);
OutputFile = fullfile(DirOutput,['M0_im_proc_' NameOutput '.jpg']);
fprintf('%s\n',['Writing ' OutputFile]);
xASL_vis_Imwrite(VisualQC, OutputFile);

end
