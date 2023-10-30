function [x] = xASL_im_CreateGroupAnalysisMask(x, Threshold)
%xASL_im_CreateGroupAnalysisMask Creates population analysis mask
%
% FORMAT: [x] = xASL_im_CreateGroupAnalysisMask(x, Threshold)
%
% INPUT:
%   x           - struct containing pipeline environment parameters, useful
%                 when only initializing ExploreASL/debugging (REQUIRED)
%   Threshold   - percentage of population that needs to have data within a
%                 voxel (OPTIONAL, DEFAULT=0.95)
%
% OUTPUT:
%   x                       - same as input
%   x.S.VBAmask             - final analysis mask, to be used later in the pipeline
%   x.S.MaskSusceptibility  - susceptibility mask, to be used later in the pipeline
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function takes the mean population-based probability
%              maps of masks, thresholds and combines them:
%
%              A. Creation GM, WM & WholeBrain masks by p>0.5
%              B. Create, combine & save vascular, susceptibity & FoV
%                 masks:
%                 - MaskVascular
%                 - MaskSusceptibility = MaskSusceptibility & MaskFoV
%              C. Create & save VBA mask
%                 - MaskAnalysis = MaskVascular & MaskSusceptibility
%                 - x.S.VBAmask = MaskAnalysis & GMmask
%              D. Visualization: Creates a figure with columns being
%                 following maps/masks overlaid over mean population T1w:
%                 1. FoV probability 0-50% missing voxels
%                 2. Vascular 0-7.5% missing voxels
%                 3. Susceptibility 0-50% missing voxels
%                 4. Analysis mask
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        xASL_im_CreateGroupAnalysisMask(x);
% __________________________________
% Copyright (c) 2015-2023 ExploreASL


%% Admin
if nargin<2 || isempty(Threshold)
    Threshold = 0.95; % default threshold
end



bSkipStandard = false;

if x.dataset.nSubjectsSessions<16
    % With too small datasets, created templated are not reliable
    fprintf('\n\n%s\n\n', ['Only n=' num2str(x.dataset.nSubjectsSessions) ' subject*runs, population-based analysis mask may not be useful']);
    % x.S.MaskSusceptibility = xASL_im_IM2Column(ones(121,145,121),x.S.masks.WBmask);
    % x.S.VBAmask = x.S.MaskSusceptibility;
end

PathTemplateSusceptibilityMask = xASL_adm_GetFileList(x.D.TemplatesStudyDir, ['^MaskSusceptibility' x.S.TemplateNumberName '_bs-mean\.nii$'], 'FPList');
PathFoV = xASL_adm_GetFileList(x.D.TemplatesStudyDir, ['^FoV' x.S.TemplateNumberName '_bs-mean\.nii$'], 'FPList', [1 1]);
PathVascularMask = xASL_adm_GetFileList(x.D.TemplatesStudyDir, ['^MaskVascular' x.S.TemplateNumberName '_bs-mean\.nii$'], 'FPList', [1 1]);

PathT1 = xASL_adm_GetFileList(x.D.TemplatesStudyDir, ['^T1' x.S.TemplateNumberName '_bs-mean\.nii$'], 'FPList', [1 1]);
PathpGM = xASL_adm_GetFileList(x.D.TemplatesStudyDir, ['^pGM' x.S.TemplateNumberName '_bs-mean\.nii$'], 'FPList', [1 1]);
PathpWM = xASL_adm_GetFileList(x.D.TemplatesStudyDir, ['^pWM' x.S.TemplateNumberName '_bs-mean\.nii$'], 'FPList', [1 1]);
PathpCSF = xASL_adm_GetFileList(x.D.TemplatesStudyDir, ['^pCSF' x.S.TemplateNumberName '_bs-mean\.nii$'], 'FPList', [1 1]);

pathsTemplates = {'PathTemplateSusceptibilityMask' 'PathFoV' 'PathVascularMask' 'PathpGM' 'PathpWM' 'PathpCSF' 'PathT1'};

for iTemplate=1:length(pathsTemplates)
    nameTemplate = pathsTemplates{iTemplate};
    pathTemplate = eval(nameTemplate);
    if isempty(pathTemplate)
        warning(['Skipping because of missing template: ' nameTemplate]);
        bSkipStandard = true;
    end
end

% Not sure if this is the same for the NativeSpaceAnalysis, this still has to be fixed by Jan.
if bSkipStandard && ~x.modules.population.bNativeSpaceAnalysis
    return;
end


if ~bSkipStandard
    % PM: these fallbacks for anatomical templates could be removed,
    % because the ASL modules won't run without anatomical images.
    % If a dataset doesn't have T1w images, these can be supplied using 
    % x.modules.asl.bUseMNIasDummyStructural == 1
    
    if isempty(PathpGM)
        warning('Missing pGM image, using MNI version');
        PathpGM = fullfile(x.D.MapsSPMmodifiedDir, 'rc1T1_ASL_res.nii');
    end
    if isempty(PathpWM)
        warning('Missing pWM image, using MNI version');
        PathpWM = fullfile(x.D.MapsSPMmodifiedDir, 'rc2T1_ASL_res.nii');
    end
    
	%% Creation GM, WM & WholeBrain masks by p>0.5
	GMmask = xASL_io_Nifti2Im(PathpGM)>0.5;
	WMmask = xASL_io_Nifti2Im(PathpWM)>0.5;
	
	if ~isempty(PathpCSF)
		pCSF = xASL_io_Nifti2Im(PathpCSF)>0.5;
		WholeBrain = (GMmask+WMmask+pCSF)>0.5;
	else
		WholeBrain = (GMmask+WMmask)>0.5;
	end

	GMmask = GMmask>0.5;
	WMmask = WMmask>0.5;
	
	%% B) Create, combine & save vascular, susceptibity & FoV masks
	MaskVascular = xASL_io_Nifti2Im(PathVascularMask)>=Threshold;
    MaskFoV = xASL_io_Nifti2Im(PathFoV)>=Threshold;
    


    % Legacy Susceptibility Masking
    MaskSusceptibility = xASL_io_Nifti2Im(PathTemplateSusceptibilityMask);
    susceptibilitySortedIntensities = sort(MaskSusceptibility(isfinite(MaskSusceptibility)));
    indexIs = round(Threshold.*length(susceptibilitySortedIntensities));
    ThresholdSuscept = susceptibilitySortedIntensities(indexIs);
    MaskSusceptibility = MaskSusceptibility >= ThresholdSuscept;

    % Combine susceptibility & FoV
    MaskSusceptibility = MaskSusceptibility & MaskFoV;    

	% Save them
	xASL_io_SaveNifti(PathFoV, fullfile(x.S.StatsDir,'MaskVascular.nii'), MaskVascular, 8, true);
	xASL_io_SaveNifti(PathFoV, fullfile(x.S.StatsDir,'MaskSusceptibility.nii'), MaskSusceptibility, 8, true);
	
	% this is used in stats:
	x.S.MaskSusceptibility = xASL_im_IM2Column(MaskSusceptibility,x.S.masks.WBmask);
end

% Define atlas paths
pathMaskInputList = ...
	{fullfile(x.S.StatsDir,'MaskSusceptibility.nii')...
	 fullfile(x.D.MapsSPMmodifiedDir,'LeftRight.nii')};
 
for iAtlas = 1:length(x.S.Atlases)	
	if strcmp(x.P.Atlas.(x.S.Atlases{iAtlas})((end-2):end),'.gz')
		pathMaskInputList{iAtlas+2} = x.P.Atlas.(x.S.Atlases{iAtlas})(1:(end-3));
	else
		pathMaskInputList{iAtlas+2} = x.P.Atlas.(x.S.Atlases{iAtlas});
	end
end

%% B2) Save FOV mask for each subject
if x.modules.population.bNativeSpaceAnalysis
	for iSession=1:x.dataset.nSessions
		%x.SESSIONS{iSession}

		% Searching for available images
		for iSubject = 1:x.dataset.nSubjects
			SubjSess = (iSubject-1)*x.dataset.nSessions + iSession;

			x.dir.SUBJECTDIR = fullfile(x.D.ROOT,x.SUBJECTS{iSubject});
			x.dir.SESSIONDIR = fullfile(x.D.ROOT,x.SUBJECTS{iSubject},x.SESSIONS{iSession});
			x = xASL_init_FileSystem(x);

			xASL_TrackProgress(SubjSess,x.dataset.nSubjects*x.dataset.nSessions);
			if xASL_exist(x.P.Path_PWI)
				x = xASL_adm_DefineASLResolution(x);
				
				for iAtlas = 1:length(pathMaskInputList)
					% Path to the atlas
					pathMaskInput = pathMaskInputList{iAtlas};
					if xASL_exist(pathMaskInput)
						% Derive the name of the mask
						[~, filenameAtlas, extensionAtlas] = xASL_fileparts(pathMaskInputList{iAtlas});
						pathMaskOutput = fullfile(x.dir.SESSIONDIR, [filenameAtlas '_Atlas' extensionAtlas]);
						
						% Automatically detect mask type
						% PM:   in the future we can also do option 1 for multi-label
						%       masks, by splitting them in multiple individual masks and
						%       treating those separately
						imMaskTmp = xASL_io_Nifti2Im(pathMaskInput);
						% Binary file
						if max(imMaskTmp(:)) == 1
							% binary masks - presmooth, spline-interpolation, cut at 50%
							% Pre-smooth the mask before downsampling to native space
							[tmpPath,tmpFile,tmpExt] = xASL_fileparts(pathMaskOutput);
							pathTmpPreSmooth = fullfile(tmpPath,['pres_' tmpFile tmpExt]);
							pathTmpPreSmooth = xASL_im_PreSmooth(x.P.Path_PWI, pathMaskInput, pathTmpPreSmooth, x.S.optimFWHM_mm,[], x.P.Path_mean_PWI_Clipped_sn_mat, 1);
							
							% Transform the Mask to native space
							xASL_spm_deformations(x, pathTmpPreSmooth, pathMaskOutput, 2, x.P.Path_PWI, x.P.Path_mean_PWI_Clipped_sn_mat, x.P.Path_y_ASL);
							xASL_delete(pathTmpPreSmooth);
							
							% Threshold at 50%
							imMaskTmp = xASL_io_Nifti2Im(pathMaskOutput);
							imMaskTmp = imMaskTmp > 0.5;
							xASL_io_SaveNifti(pathMaskOutput, pathMaskOutput, imMaskTmp);
						else % multilabel file
							% multi-label masks - no presmooth, nearest-neighbor interpolation, no thresholding
							xASL_spm_deformations(x, pathMaskInput, pathMaskOutput, 0, x.P.Path_PWI, x.P.Path_mean_PWI_Clipped_sn_mat, x.P.Path_y_ASL);
							xASL_adm_GzipNifti(pathMaskOutput);
						end
					end
				end
			end
		end
	end
	% Zip all atlases at the end
	for iAtlas = 1:length(pathMaskInputList)
		if xASL_exist(pathMaskInputList{iAtlas})
			xASL_adm_GzipNifti(pathMaskInputList{iAtlas});
		end
	end
end

if bSkipStandard
	return;
end

%% C) Create & save VBA mask
MaskAnalysis = MaskSusceptibility;
x.S.VBAmask = MaskAnalysis & GMmask; % this should be an image matrix, not a column
xASL_io_SaveNifti(PathFoV, fullfile(x.D.PopDir,'VBA_mask_final.nii'), x.S.VBAmask, 8, true);

%% D) Visualization
% 0) Admin
IntScale = [0.75 0.65];
ColorMap = {x.S.gray x.S.jet256};
bClip = [];
bWhite = false;
x.S.TraSlices = [53 35]; % 35]; % Z slice
x.S.CorSlices = 77; %
x.S.SagSlices = 63; % left right
x.S.ConcatSliceDims = 0;
x.S.bCrop = -15;
MasksAre = {ones(121,145,121) WholeBrain};


% 1) FoV probability
FoVim = 1-xASL_io_Nifti2Im(PathFoV);
FoVim(FoVim>0.5) = 0.5; % show 0-50% missing voxels in colors
ImOut{1} = xASL_vis_CreateVisualFig(x, {PathT1 FoVim}, [], IntScale, [], ColorMap, bClip, [], bWhite);

% % 2) Vascular
% ImVascular = 1-xASL_io_Nifti2Im(PathVascularMask);
% ImVascular(ImVascular>0.075) = 0.075; % show
% ImOut{2} = xASL_vis_CreateVisualFig(x, {PathT1 ImVascular}, [], IntScale, [], ColorMap, bClip, MasksAre, bWhite);

% 2) Susceptibility
ImSusceptibility = 1-xASL_io_Nifti2Im(PathTemplateSusceptibilityMask);
% ImSusceptibility(ImSusceptibility>0.5) = 0.5; % ThresholdSuscept
ImSusceptibility(ImSusceptibility<(1-ThresholdSuscept)) = 1-ThresholdSuscept;
ImSusceptibility = ImSusceptibility-(1-ThresholdSuscept);


ImOut{2} = xASL_vis_CreateVisualFig(x, {PathT1 ImSusceptibility}, [], IntScale, [], ColorMap, bClip, MasksAre, bWhite);

% 3) Resulting analysis mask
ImOut{3} = xASL_vis_CreateVisualFig(x, {PathT1 MaskAnalysis}, [], IntScale, [], ColorMap, bClip, MasksAre, bWhite);

% 5) Save image
xASL_vis_Imwrite([ImOut{1},ImOut{2},ImOut{3}], fullfile(x.S.StatsDir,'AnalysisMaskCreation.jpg'),[], false);

%     figure(1);  imshow([ImOut{1},ImOut{2},ImOut{3}],'border','tight','InitialMagnification',400)

%% Reset visualization settings
x.S.TraSlices = x.S.slicesLarge;
x.S.CorSlices = [];
x.S.SagSlices = [];

% %% Create the colorbars
% MaxValue = [0.5 0.75];
% ColorMap = ColorMap{2};
% for iBar=1:length(MaxValue)
%     if isfinite(MaxValue(iBar))
%         DummyIm = repmat([0:0.01:1].*MaxValue(iBar),[101 1]);
%         figure(iBar); imshow(DummyIm,[],'colormap',ColorMap,'InitialMagnification',400);
%         colorbar;
%     end
% end


%     Unused T1 mask
%     T1mask = xASL_io_Nifti2Im(PathT1);
%     SortInt = sort(T1mask(:));
%     Thr = SortInt(round(0.5*length(SortInt)));
%     T1mask = T1mask>Thr;


end