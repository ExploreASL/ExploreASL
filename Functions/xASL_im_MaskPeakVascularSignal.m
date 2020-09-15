function [MaskIM, TreatedCBF] = xASL_im_MaskPeakVascularSignal(PathCBF, Path_M0, bClip, ClipThresholdValue, CompressionRate)
%xASL_im_MaskPeakVascularSignal Detect voxels with high intra-vascular peak
%ASL signal
%
% FORMAT: [MaskIM, CBF] = xASL_quant_VascularContrast(PathCBF, Path_M0, CompressionRate, ClipThresholdValue, bClip)
%
% INPUT:
%   PathCBF             - path to CBF image (REQUIRED)
%   Path_M0             - path to M0 image (OPTIONAL)
%   bClip               - true for applying compression/clipping of the
%                         intra-vascular detected voxels (OPTIONAL, DEFAULT=false)
%   ClipThresholdValue  - how many Mean Absolute Difference (MAD) above the
%                         median we define as clipping/compression threshold (OPTIONAL, DEFAULT=3)
%   CompressionRate     - Compression slope of signal above
%                         ClipThresholdValue, or how "fast" the high values are compressed. Should be a
%                         value between 0 (hard clipping), 0.5 (taking root) to 1 (not affected) (OPTIONAL, DEFAULT=0.675)
%
% OUTPUT:
%   MaskIM              - Binary mask with 1 for negative signal
%   TreatedCBF          - Original CBF image, with extremely high vascular
%                         signal compressed
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function searches for an acceptable high
%              threshold as definition of high intra-vascular ASL signal.
%              It also allows to compress the values here (when
%              bClip==true). Compression retains some variability, but limits their outlying influence
%              on statistics.
%              The procedure works as follows:
%
%              1. Segment {{GM}} on {{ASL}} image at 80th percentile of {{CBF}} image distribution
%              2. For PWI & CBF images, select voxels higher than median + ClipThresholdValue {{MAD}}
%                 Vascular artifacts will have a high intensity in both images, whereas errors by division by M0 will only have a high
%                 intensity on the M0 image, and high values due to a biasfield will only be noticeable on the PWI image
%              3. Combine the two created masks
%              4. Obtain average clipping value from selected voxels from the combined masks
%              5. Apply compression if requested. If not, output image will
%              have NaNs for intra-vascular voxels.
%
%              Note that the definition of the threshold is obtained within
%              the GM only, but that this threshold is applied to the full
%              image to also remove extracranial extreme values.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: MaskIM = xASL_im_MaskPeakVascularSignal('/analysis/Sub-001/ASL_1/CBF.nii','/analysis/Sub-001/ASL_1/rM0.nii');
% __________________________________
% Copyright 2015-2020 ExploreASL




%% ===================================================================
% Admin

if nargin<1 || isempty(PathCBF)
    error('No CBF image found, skipping...');
else
    TreatedCBF = xASL_io_Nifti2Im(PathCBF);
    if sum(isfinite(TreatedCBF(:)))==0
        warning('Empty CBF image, skipping...');
        MaskIM = ones(size(TreatedCBF));
        TreatedCBF = TreatedCBF;
        return;
    end
end

if nargin<2 || isempty(Path_M0)
	M0_im = 1;
	Path_M0 = [];
else
	if ~xASL_exist(Path_M0,'file')
		warning('No usable M0 image, using CBF only');
		M0_im = 1;
	else
		M0_im = xASL_io_Nifti2Im(Path_M0);
		if ~(min(size(TreatedCBF)==size(M0_im))) % if CBF & M0 images arent the same size (which is why rM0 should be the input path in case of different sizes)
			warning('M0 & CBF images differed in size, using CBF image only');
			M0_im = 1;
		elseif xASL_stat_SumNan(M0_im(:))==0
			warning('Empty M0 image, Using CBF image only');
			M0_im = 1;
		end
	end
end

if nargin<3 || isempty(bClip)
    bClip = false; % default is set to NaN
end
if nargin<4 || isempty(ClipThresholdValue)
    ClipThresholdValue = 3; % emperically this seems a good value
end
if nargin<5 || isempty(CompressionRate)
    CompressionRate = 0.675; % emperically this seems a good value
end


%% 1) Segment GM on ASL image at 80th percentile of CBF image distribution
DataSort = sort(TreatedCBF(isfinite(TreatedCBF)));
IndexN = round(0.8*length(DataSort));
GMmask = TreatedCBF>DataSort(IndexN) & isfinite(TreatedCBF);

%% 2) For PWI & CBF images, select voxels higher than median + ClipThresholdValue MAD
% Create PWI & PWI/M0 image
TempCBF{1} = TreatedCBF.*M0_im; % PWI
TempCBF{2} = TreatedCBF; % CBF (PWI/M0)

% Get masks for both PWI & CBF images (i.e. with & without M0 division
for ii=1:2
    Mask2Use = GMmask & isfinite(TempCBF{ii});
    medianValue(ii) = median(TempCBF{ii}(Mask2Use));
    madValue(ii) = xASL_stat_MadNan(TempCBF{ii}(Mask2Use),1);
    ClipThr(ii) = medianValue(ii) + (ClipThresholdValue*madValue(ii));
    MaskIM{ii} = TempCBF{ii}>ClipThr(ii);
end

%% 3) Combine the two created masks
MaskIM = MaskIM{1} & MaskIM{2};

%% 4) Obtain average clipping value from selected voxels from the combined masks
% We compress after M0 correction, otherwise we could still have very high values
ClipThr = xASL_stat_MedianNan(TreatedCBF(MaskIM));
MaskIM = TreatedCBF>ClipThr; % create new Mask
% This will also reduce extracranial M0 division artifacts

if numel(M0_im)>1 && min(size(M0_im)==size(MaskIM)) % if they have the same image size
    MaskIM(M0_im==0 | isnan(M0_im)) = 1;
end


%% 5) Apply compression if requested
if bClip
    % Compress all values above ClipThr
    % Here, the regions with values higher than ClipThr, are compressed
    % "How much higher the values are than ClipThr is reduced by factor
    % "CompressionRate" ". We don't allow the CBF value to increase
    TreatedCBF(MaskIM) = min( TreatedCBF(MaskIM), ClipThr + ((TreatedCBF(MaskIM)-ClipThr).^CompressionRate));
else
    TreatedCBF(MaskIM) = NaN;
end


end
