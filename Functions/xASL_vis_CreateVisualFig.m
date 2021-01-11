function [ImOut, FileName] = xASL_vis_CreateVisualFig(x, ImIn, DirOut, IntScale, NamePrefix, ColorMap, bClip, MaskIn, bWhite, MaxWindow, bTransparancy, bVerbose)
% xASL_vis_CreateVisualFig Flexible tool to create figure for visualization
% of standard space images
%
% FORMAT: [ImOut, FileName] = xASL_vis_CreateVisualFig(x, ImIn, DirOut, IntScale, NamePrefix, ColorMap, bClip)
%
% INPUT:
%   x            - structure containing fields with information when this function is called from ExploreASL toolbox (OPTIONAL)
%   ImIn         - input image, can be a path to a NIfTI file or an image matrix. Can be a single image, or multiple images,
%                  in the latter case the other images are overlaid over the first (e.g. to visualize segmentation).
%                  Input images should be in [1.5 1.5 1.5] mm MNI space
%                  The same number of inputs should be provided to IntScale and ColorMap (REQUIRED)
%   DirOut       - folder where output Figure file should be written (OPTIONAL, DEFAULT=omit writing output Figure file)
%   IntScale     - vector specifying the relative intensity of each image layer, with same length of as ImIn (OPTIONAL, DEFAULT=1)
%                  this allows setting transparancy, when the intensity of the overlay is reduced
%   NamePrefix   - prefix for filename (OPTIONAL, DEFAULT = '')
%   ColorMap     - colormaps used for each layer (OPTIONAL, DEFAULT =x.S.colors_ROI
%                  x.S.colors_ROI = {gray red blue green yellow purple turqoise orange}
%   bClip        - vector, false for disabling clipping (OPTIONAL, DEFAULT=true for first image, false for overlays)
%   MaskIn       - cell structure containing 2 binary masks matrices (OPTIONAL, DEFAULT=no masking)
%   bWhite       - true for switching background to white (OPTIONAL, DEFAULT=black background)
%   MaxWindow    - cell structure containing 2 maximal values (ceiling clipping) (OPTIONAL, DEFAULT=automatic window-leveling)
%   bTransparancy - true for transparant results when overlaying a mask (OPTIONAL, DEFAULT=false)
%   bVerbose      - true for feedback on what this function does (OPTIONAL, DEFAULT=false)
%
% INPUT FIELDS IN x, USED BY xASL_vis_TransformData2View:
%              ORIENTATION SETTINGS:
%               x.S.TraSlices - which transversal slices to show (n, orientation omitted if empty, DEFAULT = 20:7:97)
%               x.S.CorSlices - which coronal slices to show (n, orientation omitted if empty, DEFAULT = empty)
%               x.S.SagSlices - which sagittal slices to show (n, orientation omitted if empty, DEFAULT = empty)
%              CROP SETTINGS:
%               x.S.bCrop - number of voxels you want to crop from the default size, single value,
%                           ranging from 50 (pieces of brain cut out) to -Inf (creating space between images) (DEFAULT=0)
%              CONCATENATION SETTINGS:
%              x.S.ConcatSliceDims % whether orientations should be concatenated horizontally (1, DEFAULT) or vertically (0)
%              x.S.Square % whether we prefer to tile square (DEFAULT = true)
%
% OUTPUT:
%   ImOut        - output image matrix, in Matlab colorscale (3 images for RGB) (OPTIONAL)
% OUTPUT FILE    - saved figure (e.g. as jpg), same as ImOut
%
% DESCRIPTION: This function creates a visualization Figure by merging flexibly rearranging NIfTI slices, input matrix or
%              path, managing colormaps for different merged image layers. Current use is for visual QC figures and overview in papers.
%              Function is structured as:
%
%              1. Admin, deal with input arguments
%              2. Process image layers separately
%                 * xASL_im_TransformData2View: Reshapes image data into visualization figure
%                 * xASL_im_ClipExtremes: Clips image to given percentile
%                   also we scale for peak intensity, we make sure that there is no
%                   visible clipping/distortion
%                 * Convert to colors, using any input colormaps
%              3. combine image layers, using input argument IntScale
%              4. print figure
%
%              This function assumes that the first image is a grayscale background
%              image (e.g. for transparancy reasons), if there are multiple
%              images
%
% EXAMPLE: for overlaying red GM segmentation over grayscale T1w background:
%          [ImOut] = xASL_vis_CreateVisualFig(x, {'//AnalysisDir/Population/T1_Sub-001.nii' '//AnalysisDir/Population/rc1T1_Sub-001.nii', [], [1 0.5], [], {x.S.gray x.S.red});
%          or for printing to file:
%          xASL_vis_CreateVisualFig(x, {'//AnalysisDir/Population/T1_Sub-001.nii' '//AnalysisDir/Population/rc1T1_Sub-001.nii', '//AnalysisDir/Population/T1w_Check', [1 0.5], 'pGM', {x.S.gray x.S.red});
% __________________________________
% Copyright 2015-2020 ExploreASL



%% ----------------------------------------------
%% 1. Admin
ImOut = NaN; % if anything goes wrong

% Get visualization settings if required
if nargin<1 || isempty(x)
    x = xASL_init_VisualizationSettings;
elseif ~isfield(x,'S') || isempty(x.S) || ~isfield(x.S, 'TransCrop') || ~isfield(x.S, 'jet256')
    x = xASL_init_VisualizationSettings(x);
end
if ~iscell(ImIn) % if single input image/path, this can be inside/outside cell
    ImIn = {ImIn};
end
if nargin<3 || isempty(DirOut)
    DirOut=NaN;
end
if nargin<4 || isempty(IntScale)
    IntScale = ones(1,length(ImIn));
end
if nargin<5 || isempty(NamePrefix)
    NamePrefix = '';
end
if nargin<6 || isempty(ColorMap)
    ColorMap = x.S.colors_ROI;
elseif ~iscell(ColorMap)
    ColorMap = {ColorMap};
end
if nargin<7 || isempty(bClip)
    bClip = zeros(1,length(ImIn));
    bClip(1) = true; % by default clip first image, leave overlays unclipped
end
if length(bClip)<length(ImIn) % They need to have the same length - add zeros to the end
	bClip((length(bClip)+1):length(ImIn)) = 0;
end
if nargin<8 || isempty(MaskIn)
    for iIm = 1:length(ImIn)
        MaskIn{iIm} = ones(121,145,121); % default = no masking
    end
elseif ~iscell(MaskIn)
    MaskIn = {MaskIn};
end
if nargin<9 || isempty(bWhite)
    bWhite = false; % default = black background
end
if nargin<10 || isempty(MaxWindow)
    for iIm = 1:length(ImIn)
        MaxWindow{iIm} = 0; % default = automatic window leveling
    end
end
if nargin<11 || isempty(bTransparancy)
    bTransparancy = false;
end
if nargin<12 || isempty(bVerbose)
    bVerbose = false;
end

if isempty(ImIn)
    fprintf('%s\n', 'xASL_vis_CreateVisualFig:No image data, skipping image creation');
    FigureOut{1} = NaN;
    return;
end


% First check file existence (only for input images provided as FilePath)
for iC=1:length(ImIn)
    if  numel(ImIn{iC})<512 % assume this is a FilePath
        if ~xASL_exist(ImIn{iC}, 'file')
            fprintf('%s\n','xASL_vis_CreateVisualFig: No image data, skipping image creation');
            FigureOut{1} = NaN;
            return;
        else
            [~, Ffile{iC}] = xASL_fileparts(ImIn{iC});
        end
    end
end



%% ------------------------------------------------------------------------
%% 2. Process image layers
for iIm=1:length(ImIn)
    IM{iIm} = xASL_io_Nifti2Im(ImIn{iIm}); % load if NIfTI, pass through if matrix (IM)

    if length(MaskIn)<iIm
        MaskIn{iIm} = ones(121,145,121); % default=no masking
    end

    DimIm = size(IM{iIm});

    if length(DimIm)>3
        fprintf('Minor warning, xASL_vis_CreateVisualFig: image had multiple 3D images, using first only\n');
    end
    IM{iIm} = IM{iIm}(:,:,:,1,1,1,1,1,1,1,1,1);

    if max(size(IM{iIm})~=[121 145 121])
        % resample first
        warning('Input image didnt have correct dims, resample first to standard space for visualization, skipping...');
        FigureOut{1} = NaN;
        return;
    end

    %% 0. Masking & shifting to zero
    MaskIn{iIm}(isnan(MaskIn{iIm})) = 0; % Remove NaNs
    IM{iIm}(~logical(MaskIn{iIm})) = 0;

    %% First clipping
    if bClip(iIm)
        IM{iIm}(IM{iIm}<0) = 0;
    end
    
    % Shift to zero on histogram
    MinIntIm = min(IM{iIm}(:));
    if MinIntIm>0
        IM{iIm} = IM{iIm}-MinIntIm;
    end

    %% a. transforms the input images to a visualization format
    % first do this for the mask
    if bWhite
        % do the same for the mask, allowing to switch background to white
        if xASL_stat_SumNan(MaskIn{iIm}(:)==0)==0 % in case there is no mask
            % create a mask from the data
            CornerMask = logical(zeros(size(IM{iIm})));
            CornerMask(1,1,1) = 1;
            CornerMask(end,1,1) = 1;
            CornerMask(1,end,1) = 1;
            CornerMask(end,end,1) = 1;
            CornerMask(1,1,end) = 1;
            CornerMask(end,1,end) = 1;
            CornerMask(1,end,end) = 1;
            CornerMask(end,end,end) = 1;

            Regions = spm_bwlabel(double(IM{iIm}==0));
            if sum(Regions(:))>0
                % get the largest region
                for iR=1:length(Regions)
                    SizeRegion(iR,1)=sum(sum(sum(Regions==iR)));
                end
                IndexLargestRegion = find(SizeRegion==max(SizeRegion(:)));
                % if all corners lie in this region, make this the mask
                if sum(sum(sum(Regions==IndexLargestRegion & CornerMask)))==8
                    MaskIn{iIm} = Regions~=IndexLargestRegion;
                else
                    MaskIn{iIm} = logical(ones(121,145,121));
                    % keep the mask full, we cannot mask
                end
            else
                MaskIn{iIm} = logical(ones(121,145,121));
                % keep the mask full, we cannot mask
            end
        end
        MaskWhite{iIm} = xASL_vis_TransformData2View(MaskIn{iIm}(:,:,:,1), x);
    end

    % now do this for the data
    IM{iIm}(isnan(IM{iIm})) = 0; % Remove NaNs
    IM{iIm} = xASL_vis_TransformData2View(IM{iIm}(:,:,:,1), x);

    %% b. Clipping image intensities for image layer
    if MaxWindow{iIm}~=0
        MaxInt1 = MaxWindow{iIm};
        MaxInt2 = MaxWindow{iIm};
    else
        MaxInt1 = 0.99; % default
        MaxInt2 = 0.975;
    end

    if bClip(iIm)
        IM{iIm} = xASL_im_ClipExtremes(IM{iIm}, MaxInt1, 0, 0);
        % this function clips the image, to avoid a too large viewing window
    end

    Um = unique(IM{iIm}(:));
    if length(Um)==1 && Um==0
        warning('This image cannot be sorted, skipping');
        return;
    end
        
    SortInt = sort(IM{iIm}(:));
    MaxInt2 = SortInt(round(MaxInt2*length(SortInt)));

    if bVerbose && iIm==1
        fprintf('%s\n', ['Min & max value were [' xASL_num2str(min(IM{iIm}(:))) ' ' xASL_num2str(MaxInt2) '] for background image']);
    elseif bVerbose
        fprintf('%s\n', ['Min & max value were [' xASL_num2str(min(IM{iIm}(:))) ' ' xASL_num2str(MaxInt2) '] for overlay image ' xASL_num2str(iIm-1)]);
    end

    if bClip(iIm)
        IM{iIm}(IM{iIm}>MaxInt2) = MaxInt2;
        IM{iIm} = IM{iIm} ./MaxInt2;
        % Here we set the 5th upper percentile as maximum intensity,
        % this is an arbitrary windowing level setting that often works
    end

    %% c) Convert image layer to colors
    SzC = size(ColorMap{iIm},1);
    if SzC<250 || SzC>260
        warning(['Colormap expects 256 indices, but had only ' num2str(SzC)]);
        fprintf('This can lead to image clipping artifacts\n');
        % PM: here we can resample the input colormap into 256 if it is a power of
        % e.g. 32/64
    end

    IM{iIm} = IM{iIm}./max(IM{iIm}(:)); % divide image by max to avoid color clipping
    IM{iIm} = ind2rgb(round(IM{iIm}.*255), ColorMap{iIm}); % convert the image to a RGB false colored-image
    IM{iIm} = IM{iIm} .* IntScale(iIm); % scale the image, if requested
end

%% 3. Combine image layers
% We can add the image layers with the intensities set by IntScale above.
% However, if we overlay a binary mask on top of a grayscale image with
% high intensities, this will go wrong: then we need to force transparancy

ImOut = IM{1};
if length(IM)>1 % if we have more images, we overlay them on the first
    TransparancyMask = sum((IM{1} ./ IntScale(1))>0.75,3)==3;
    TransparancyMask = repmat(TransparancyMask,[1 1 3]);
    % TransparancyMask contains very bright voxels, where we may need transparancy

    Mask2 = zeros(size(TransparancyMask));
    for iImage=2:length(IM)
        NonGrayMask = ~(IM{iImage}(:,:,1)==IM{iImage}(:,:,2) & IM{iImage}(:,:,2)==IM{iImage}(:,:,3));
        NonGrayMask = repmat(NonGrayMask,[1 1 3]);
        % Create mask where overlaid image has non-zero non-gray values
        NonGrayMask = IM{iImage}>0 & NonGrayMask;
        % Take all voxels within this mask that contain any RGB color
        Mask2 = Mask2 | NonGrayMask;
        % Combine this mask for all image layers
    end

    TransparancyMask = TransparancyMask & Mask2; % combine masks

    IsBinaryMask = true;
    if ~length(unique(IM{iImage}(:)))==2 % if it is a binary mask
        IsBinaryMask = false;
    end
    if IsBinaryMask && bTransparancy
        IntFactor = IntScale./sum(IntScale); % do transparancy
    else
        IntFactor(length(IntScale)) = 1;
        IntFactor(1) = 1;
    end

    ImOut(TransparancyMask) = ImOut(TransparancyMask).*IntFactor(1);

    for iImage=2:length(IM)
        ImOut(TransparancyMask) = ImOut(TransparancyMask) + IM{iImage}(TransparancyMask).*IntFactor(iImage);
        ImOut(~TransparancyMask) = ImOut(~TransparancyMask) + IM{iImage}(~TransparancyMask);
    end

    %% PM: create a transparancy setting
end

%% 3.5. Deal with background color
if bWhite
    % first create common mask, by taking the MIP, a voxel will only be
    % masked out (i.e. set to white), if it doesnt exist in any of the
    % masks
    for iMask=1:length(MaskWhite)
        if iMask==1
            CombinedMask = MaskWhite{iMask};
        else
            CombinedMask = max(CombinedMask, MaskWhite{iMask}); % empty masks are ignored
        end
    end
    % convert mask to "color mask"
    CombinedMask = repmat(logical(CombinedMask),[1 1 3]);
    ImOut = ImOut./max(ImOut(:)); % rescaling needed to set background to white instead of gray
    ImOut(~CombinedMask) = 1; % set background to white
end

%%  -----------------------------------------------------------------
%%  4. Print Figure
% the output image is always passed through as output argument, if
% requested. If DirOut is provided, this image is also saved as jpg file
if numel(DirOut)>1 && numel(ImOut)>1 % when DirOut or ImOut==NaN, skip this
    xASL_adm_CreateDir(DirOut);
    for iC=1:length(Ffile)
        FileName = [NamePrefix Ffile{iC}];
    end
    OutputFile = fullfile(DirOut,[FileName '.jpg']);
    xASL_vis_Imwrite(ImOut, OutputFile);
end

if ~exist('FigureOut', 'var')
    FigureOut{1} = NaN;
end


end
