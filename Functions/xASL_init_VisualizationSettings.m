function [x] = xASL_init_VisualizationSettings(x)
%xASL_init_VisualizationSettings Defines the visualization settings for ExploreASL
%
% FORMAT:   [x] = xASL_init_VisualizationSettings(x)
%
% INPUT:
%  x 	    - structure containing fields with information when this function is called from within the ExploreASL pipeline (REQUIRED)
%
% OUTPUT:
%  x 	    - structure containing fields with information when this function is called from within the ExploreASL pipeline (REQUIRED)
%
% DESCRIPTION:  This function defines several visualization settings are
%               used throughout ExploreASL's pipeline and tools, assuming a [121 145 121]
%               matrix with 1.5 mm isotropic resolution in MNI space.
%
% EXAMPLE:      [x] = xASL_init_VisualizationSettings(x);
%
% INCLUDING:
%   Slices*    - defines which transversal slices to show by default
%   TransCrop  - defines default cropping settings for all orientations
%   ColorMaps  - Generate several colormaps
%   VoxelSize  - default voxel-size in standard space/MNI as used by ExploreASL: 1.5 mm
%   skull      - brainmask for skullstripping
%
%
% __________________________________
% Copyright 2015-2020 ExploreASL


    %% Admin
    if nargin<1
        warning('Input argument "x" struct missing');
        x = struct;
    end

    if ~isfield(x,'S')
        x.S = struct;
    end

    if ~isfield(x.S,'bMasking') || isempty(x.S.bMasking)
        x.S.bMasking = [1 1 1 1];
    elseif isequal(x.S.bMasking, 1)
        x.S.bMasking = [1 1 1 1];
    elseif isequal(x.S.bMasking, 0)
        x.S.bMasking = [0 0 0 0];        
    end

    %% Slice numbers
    % Defines which transversal slices to use by default
    x.S.slices           = [53 62 74 87]; % for 1.5mm MNI, ([85 102 119 131] would be the same for 1mm MNI)
    x.S.slicesLarge      = 20:7:97;
    x.S.slicesExtraLarge = 20+([1:25]-1).*round((97-20)/24);

    x.S.nSlices          = length(x.S.slices);
    x.S.nSlicesLarge     = length(x.S.slicesLarge);
    x.S.nSlicesExtraLarge= length(x.S.slicesExtraLarge);


    %% Cropping settings
    % defines default cropping settings for all orientations
    MatrixSize = [121 145 121]; % for MNI 1.5x1.5x1.5 mm
    % This is for A P L R x.S.S I
    % 10 cropped at all sizes (assumed that everything is already in standard
    % space. Superior & Inferior little bit different to accommodate viewing
    % it next to other slices.

    x.S.TransCrop      = [10 MatrixSize(2)-10 -2 MatrixSize(1)+2 -2 123]; % this works for all orientations

    %% Create colormaps
    jet256(  1: 32,:)    = [zeros(32,1) zeros(32,1) [0.5+1/64:1/64:1]'];
    jet256( 33: 96,:)    = [zeros(64,1) [1/64:1/64:1]' ones(64,1)];
    jet256( 97:160,:)    = [[1/64:1/64:1]' ones(64,1) [1-1/64:-1/64:0]'];
    jet256(161:224,:)    = [ones(64,1) [1-1/64:-1/64:0]' zeros(64,1)];
    jet256(225:256,:)    = [[1-1/64:-1/64:0.5]' zeros(32,1) zeros(32,1)];
    jet256(1,:)          = 0;

    x.S.jet256           = jet256;

    x.S.gray             = repmat([0:1/255:1]',[1 3]);
    x.S.red              = [[0:1/255:1]' zeros(256,1) zeros(256,1)];
    x.S.yellow           = [[0:1/255:1]' [0:1/255:1]' zeros(256,1)];
    x.S.green            = [zeros(256,1) [0:1/255:1]' zeros(256,1)];
    x.S.blue             = [zeros(256,1) zeros(256,1) [0:1/255:1]'];
    x.S.purple           = [[0:1/255:1]' zeros(256,1) [0:1/255:1]'];
    x.S.turqoise         = [zeros(256,1) [0:1/255:1]' [0:1/255:1]'];
    x.S.orange           = [[0.5+1/512:1/512:1]'  [0:1/255:1]' zeros(256,1)];

    x.S.colors_ROI{1} = x.S.gray;
    x.S.colors_ROI{2} = x.S.red;
    x.S.colors_ROI{3} = x.S.blue;
    x.S.colors_ROI{4} = x.S.green;
    x.S.colors_ROI{5} = x.S.yellow;
    x.S.colors_ROI{6} = x.S.purple;
    x.S.colors_ROI{7} = x.S.turqoise;
    x.S.colors_ROI{8} = x.S.orange;

    % Generate cool, colorbar decrease
    x.S.cool = [[1/256:(1/256):1]' [1-1/256:-1/256:0]' ones(256,1)];

    % Generate hot, colorbar increase
    x.S.hot(  1: 96,:) = [[1/96:(1/96):1]' zeros(96,1) zeros(96,1)];
    x.S.hot( 97:192,:) = [ones(96,1) [1/96:(1/96):1]'  zeros(96,1)];
    x.S.hot(193:256,:) = [ones(64,1) ones(64,1)        [1/64:(1/64):1]'];

    % Adapt to make sure that ends of spectrum are not completely white or black, since this is difficult to show on a grayscale (e.g. mean T1) background
    x.S.hot = squeeze(xASL_im_ResampleLinear(x.S.hot(49:224,:), [256 3]));
    x.S.cool = squeeze(xASL_im_ResampleLinear(x.S.cool(49:256,:), [256 3]));

    x.S.VoxelSize = 1.5; % mm voxel-size of reslicing & DARTEL. DEFAULT = 1.5 mm


    %% Loading MNI brainmask
    % tight (dichotomous) mask
    if ~isfield(x,'D')
        x.D = struct;
    end
    if ~isfield(x.D,'TemplateDir') || isempty(x.D.TemplateDir)
        warning('No x.D.TemplateDir field found, not defining brainmask');
    else
        x.S.masks.skull = logical(xASL_io_Nifti2Im(fullfile(x.D.MapsSPMmodifiedDir, 'brainmask.nii')));
    end

    if ~isfield(x,'S') % Make dummy variable "S", if it doesn't exist
        x.S.SubjectSessionID = '';
        x.S.SetsID = [];
        x.S.SetsName = '';
        x.S.SetsOptions = '';
        x.S.Sets1_2Sample = [];
    end

    if ~isfield(x,'BILAT_FILTER')
        % The artifact for which this filter is intended doesn't
        % occur frequently, so don't produce an error when user has not
        % specified whether or not to run this filter.
        x.settings.BILAT_FILTER = false;
    end

    ImageWB = xASL_io_Nifti2Im(fullfile(x.D.MapsSPMmodifiedDir,'WholeBrain.nii'));
    if x.S.bMasking(4)==0
        x.S.masks.WBmask = true(size(ImageWB));
    else
        x.S.masks.WBmask = logical(ImageWB);
    end

    % Default is high quality
    if ~isfield(x.settings,'Quality')
        x.settings.Quality = 1; 
    end

    warning('off','images:initSize:adjustingMag'); % warning about scaling if an image doesnt fit screen, disable


    %% Create 64 label colors
    LabelPath = fullfile(x.D.MapsSPMmodifiedDir, 'LabelColors.mat');
    if ~exist(LabelPath, 'file')
        xASL_init_CreateLabelColors(LabelPath, x);
    else
         LabelClr = load(LabelPath,'LabelClr');
         x.S.LabelClr = LabelClr.LabelClr;
    end

end

function [x] = xASL_init_CreateLabelColors(LabelPath, x)
%xASL_init_VisualizationSettings Defines the visualization settings for ExploreASL

    % Create predefined label colors
    LabelClr(1,:)   = [0   1 1]; % Turqoise
    LabelClr(2,:)   = [1   0 0]; % Red
    LabelClr(3,:)   = [1   1 0]; % Yellow
    LabelClr(4,:)   = [0.5   0.5 1];  % Blue, slightly lighter, otherwise difficult to spot
    LabelClr(5,:)   = [0   1 0]; % Green
    LabelClr(6,:)   = [1   0 1]; % Purple
    LabelClr(7,:)   = [1 1/3 0]; % Orange
    LabelClr(8,:)   = [2/3 1 1]; % LightGreen
    LabelClr        = round(LabelClr.*5)./5; % avoid creating colors that are similar

    % Create list of colors
    ColorNfactor = 5;
    ColorsFactor = [0:1/ColorNfactor:1];
    NextN = 1;
    for i1=ColorsFactor
        for i2=ColorsFactor
            for i3=ColorsFactor
                SortC(NextN,:) = [i1 i2 i3];
                NextN = NextN+1;
            end
        end
    end

    % Pseudo-randomizing resort this list
    SortC(:,4) = randn(size(SortC,1),1)+100;
    SortC = sortrows(SortC,4);
    SortC = SortC(:,1:3);

    % keep on adding colors to LabelClr
    for iC=1:size(SortC,1)
        % Check all existing colors whether the new added one doesn't
        % already exist
        ColorExist = 0;
        for iL=1:size(LabelClr,1)
            if  min(LabelClr(iL,:)==SortC(iC,:))
                ColorExist = 1; % skip this random color, since it already exists
            end
        end

        if ~ColorExist
            LabelClr(end+1,:) = SortC(iC,:);
        end
    end

    save(LabelPath, 'LabelClr');
    % PM: Create more options by randomly permuting options, rounding them into
    % a meaningful color, checking whether it hasn't been used before &
    % creating those as separate label colors

end
