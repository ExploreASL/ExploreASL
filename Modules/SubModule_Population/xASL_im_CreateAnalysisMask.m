function [x] = xASL_im_CreateAnalysisMask(x, Threshold)
%xASL_im_CreateAnalysisMask Creates population analysis mask
%
% FORMAT: [x] = xASL_im_CreateAnalysisMask(x, Threshold)
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
%   x.S.MaskVascular        - vascular mask, to be used later in the pipeline
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
% EXAMPLE: xASL_im_CreateAnalysisMask(x);
% __________________________________
% Copyright 2015-2020 ExploreASL


%% Admin
if nargin<2 || isempty(Threshold)
    Threshold = 0.95; % default threshold
end

PathSusceptibilityMask = xASL_adm_GetFileList(x.D.TemplatesStudyDir, ['^MaskSusceptibility_n' xASL_num2str(x.nSubjectsSessions) '_bs-mean\.nii$'], 'FPList');

bSkipStandard = 0;

if x.nSubjectsSessions<16
    fprintf('%s\n', ['Too few subjects (' num2str(x.nSubjectsSessions) ') to create population-based analysis mask']);
    x.S.MaskSusceptibility = xASL_im_IM2Column(ones(121,145,121),x.WBmask);
    x.S.VBAmask = x.S.MaskSusceptibility;
    bSkipStandard = 1;
elseif isempty(PathSusceptibilityMask) && ~strcmpi(x.Sequence,'3d_spiral')
    warning('Missing susceptibility maps, skipping...');
    bSkipStandard = 1;
end

if bSkipStandard && ~x.bNativeSpaceAnalysis
    return;
    % we skip creating a group mask if we have fewer images than specified
    % above.
    % Not sure if this is the same for the NativeSpaceAnalysis, this still
    % has to be fixed by Jan
end

% Define paths to create
VBA_maskPath = fullfile(x.D.PopDir,'VBA_mask_final.nii'); % later move to statsdir
MaskVascularPath = fullfile(x.S.StatsDir,'MaskVascular.nii');
MaskSusceptibilityPath = fullfile(x.S.StatsDir,'MaskSusceptibility.nii');

% Define pre-existing paths, including warning when less or more than 1 are found
% First for SubjectsSessions (e.g. ASL)
PathFoV = xASL_adm_GetFileList(x.D.TemplatesStudyDir, ['^FoV_n' xASL_num2str(x.nSubjectsSessions) '_bs-mean\.nii$'], 'FPList', [1 1]);
PathVascularMask = xASL_adm_GetFileList(x.D.TemplatesStudyDir, ['^MaskVascular_n' xASL_num2str(x.nSubjectsSessions) '_bs-mean\.nii$'], 'FPList', [1 1]);
% Then for Subjects (e.g. structural)
PathpGM = xASL_adm_GetFileList(x.D.TemplatesStudyDir, ['^pGM_n' xASL_num2str(x.nSubjects) '_bs-mean\.nii$'], 'FPList', [1 1]);
PathpWM = xASL_adm_GetFileList(x.D.TemplatesStudyDir, ['^pWM_n' xASL_num2str(x.nSubjects) '_bs-mean\.nii$'], 'FPList', [1 1]);
PathpCSF = xASL_adm_GetFileList(x.D.TemplatesStudyDir, ['^pCSF_n' xASL_num2str(x.nSubjects) '_bs-mean\.nii$'], 'FPList', [1 1]);
PathT1 = xASL_adm_GetFileList(x.D.TemplatesStudyDir, ['^T1_n' xASL_num2str(x.nSubjects) '_bs-mean\.nii$'], 'FPList', [1 1]);

if ~isempty(PathSusceptibilityMask); PathSusceptibilityMask = PathSusceptibilityMask{1}; end
if ~isempty(PathFoV); PathFoV = PathFoV{1}; end
if ~isempty(PathVascularMask); PathVascularMask = PathVascularMask{1}; end
if ~isempty(PathpGM); PathpGM = PathpGM{1}; end
if ~isempty(PathpWM); PathpWM = PathpWM{1}; end
if ~isempty(PathpCSF); PathpCSF = PathpCSF{1}; end
if ~isempty(PathT1); PathT1 = PathT1{1}; end


if ~bSkipStandard
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
	
	DoSusceptibility = true;
	if ~isfield(x,'Sequence')
		error('x.Sequence parameter missing!');
	elseif strcmpi(x.Sequence,'2D_EPI')
		Path_Template = fullfile(x.D.MapsDir,'Templates','Susceptibility_pSignal_2D_EPI.nii');
	elseif strcmpi(x.Sequence,'3D_GRASE')
		Path_Template = fullfile(x.D.MapsDir,'Templates','Susceptibility_pSignal_3D_GRASE.nii');
	else
		% don't mask susceptibility artifacts
		DoSusceptibility = false;
	end

	MaskFoV = xASL_io_Nifti2Im(PathFoV)>=Threshold;
	if DoSusceptibility
        MaskSusceptibility = xASL_io_Nifti2Im(PathSusceptibilityMask);        

        TemplateMask = xASL_io_Nifti2Im(Path_Template)~=1; % MaskSusceptibility<0.85 &
		TemplateMask(:,:,1:10) = 0;
		TemplateMask = TemplateMask & MaskSusceptibility<0.8;
		ThrValues = MaskSusceptibility(TemplateMask);
		ThresholdSuscept = Threshold.*max(ThrValues(:));

		% Remove some extremes using 95% percentile,
		% as here we aim to obtain the 95% threshold within the probability map
		% (i.e. excluding the manually created mask to isolate the skull
		% regions of air cavities)
		MaskSusceptibility = MaskSusceptibility > ThresholdSuscept;

		% Combine susceptibility & FoV
		MaskSusceptibility = MaskSusceptibility & MaskFoV;
	else
		MaskSusceptibility = MaskFoV; % if we don't have susceptibility artifacts
	end

	% Save them
	xASL_io_SaveNifti(PathFoV, MaskVascularPath, MaskVascular, 8, true);
	xASL_io_SaveNifti(PathFoV, MaskSusceptibilityPath, MaskSusceptibility, 8, true);
	% this is used in stats:
	x.S.MaskSusceptibility = xASL_im_IM2Column(MaskSusceptibility,x.WBmask);
end
%% B2) Save FOV mask for each subject
if x.bNativeSpaceAnalysis
	for iSession=1:x.nSessions
		%x.SESSIONS{iSession}

		% Searching for available images
		for iSubject = 1:x.nSubjects
			SubjSess = (iSubject-1)*x.nSessions + iSession;

			x.SUBJECTDIR = fullfile(x.D.ROOT,x.SUBJECTS{iSubject});
			x.SESSIONDIR = fullfile(x.D.ROOT,x.SUBJECTS{iSubject},x.SESSIONS{iSession});
			x = xASL_init_FileSystem(x);

			xASL_TrackProgress(SubjSess,x.nSubjects*x.nSessions);
			if xASL_exist(x.P.Path_PWI)
				x = xASL_adm_DefineASLResolution(x);
				listMasks = {MaskSusceptibilityPath fullfile(x.D.MapsSPMmodifiedDir,'TotalGM.nii')...
					fullfile(x.D.MapsSPMmodifiedDir,'DeepWM.nii') fullfile(x.D.MapsSPMmodifiedDir,'MNI_Structural.nii')...
					fullfile(x.D.MapsSPMmodifiedDir,'LeftRight.nii') fullfile(x.D.AtlasDir,'Hammers.nii')...
                    fullfile(x.D.AtlasDir,'HOcort_CONN.nii') fullfile(x.D.AtlasDir,'HOsub_CONN.nii')};
				listOutputs = {x.P.Path_MaskSusceptibilityPop x.P.Path_TotalGMPop x.P.Path_DeepWMPop x.P.Path_MNIStructuralPop x.P.Path_LeftRightPop x.P.Path_HammersPop x.P.Path_HOcort_CONNPop x.P.Path_HOsub_CONNPop};
				MaskType  = [1 1 1 2 2 2 2 2];
				% 1 - binary masks - presmooth, spline-interpolation, cut at 50%
				% 2 - multi-label masks - no presmooth, nearest-neighbor interpolation, no thresholding
                % PM:   in the future we can also do option 1 for multi-label
                %       masks, by splitting them in multiple individual masks and
                %       treating those separately

				for kk = 1:length(listMasks)
					if xASL_exist(listMasks{kk},'file')
						switch (MaskType(kk))
							case 1
								% Pre-smooth the mask before downsampling to native space
								[tmpPath,tmpFile,tmpExt] = xASL_fileparts(listOutputs{kk});
								pathTmpPreSmooth = fullfile(tmpPath,['pres_' tmpFile tmpExt]);
								pathTmpPreSmooth = xASL_im_PreSmooth(x.P.Path_PWI,listMasks{kk},pathTmpPreSmooth,x.S.optimFWHM_mm,[],x.P.Path_mean_PWI_Clipped_sn_mat, 1);

								% Transform the Mask to native space
								xASL_spm_deformations(x, pathTmpPreSmooth, listOutputs{kk}, 2, x.P.Path_PWI, x.P.Path_mean_PWI_Clipped_sn_mat, x.P.Path_y_ASL);
								xASL_delete(pathTmpPreSmooth);

								% Threshold at 50%
								imMaskTmp = xASL_io_Nifti2Im(listOutputs{kk});
								imMaskTmp = imMaskTmp > 0.5;
								xASL_io_SaveNifti(listOutputs{kk},listOutputs{kk},imMaskTmp);
							case 2
								xASL_spm_deformations(x, listMasks{kk}, listOutputs{kk}, 0, x.P.Path_PWI, x.P.Path_mean_PWI_Clipped_sn_mat, x.P.Path_y_ASL);
								xASL_adm_GzipNifti(listOutputs{kk});
						end
					end
				end
			end
		end
	end
	for kk = 1:length(listMasks)
		if xASL_exist(listMasks{kk},'file')
			xASL_adm_GzipNifti(listMasks{kk});
		end
	end
end

if bSkipStandard
	return;
end

%% C) Create & save VBA mask
MaskAnalysis = MaskSusceptibility;
x.S.VBAmask = MaskAnalysis & GMmask; % this should be an image matrix, not a column
xASL_io_SaveNifti(PathFoV, VBA_maskPath, x.S.VBAmask, 8, true);

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
if DoSusceptibility
    ImSusceptibility = 1-xASL_io_Nifti2Im(PathSusceptibilityMask);
    % ImSusceptibility(ImSusceptibility>0.5) = 0.5; % ThresholdSuscept
    ImSusceptibility(ImSusceptibility<(1-ThresholdSuscept)) = 1-ThresholdSuscept;
    ImSusceptibility = ImSusceptibility-(1-ThresholdSuscept);
    ImSusceptibility(~TemplateMask) = 0;
else
    ImSusceptibility = MaskSusceptibility;
end

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