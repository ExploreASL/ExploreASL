function LesionIM = xASL_im_Lesion2Mask(LesionPath, x)
% xASL_im_Lesion2Mask Create multiple masks from single ROI/lesion mask, to
% be used as custom "atlas"
%
% FORMAT: LesionIM = xASL_im_Lesion2Mask(LesionPath, x)
% 
% INPUT:
%   LesionPath  - (string) path to the NIfTI containing the mask of
%                  ROI/lesion (REQUIRED)
%   x           - structure containing fields with all information required to run this function within ExploreASL (REQUIRED)
%
% OUTPUT:       - new image containing the masks
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function takes a mask and adds several ROIs, to be used as custom "atlas", e.g. when computing region-average CBF values.
% The mask % can be an ROI or lesion, if we assume it is a lesion, the following masks are created:
%
% 1. Intralesional
% 2. Perilesional (15 mm rim around the lesion)
% 3. Hemisphere (ipsilateral to lesion)
% 4. Contralateral version of 1
% 5. Contralateral version of 2
% 6. Contralateral version of 3
%
% All these masks are masked by a brainmask (pGM+pWM)>0.5
% 
% This function performs the following steps:
%
% 1. If lesion is empty, skip this & delete the file
% 2. BrainMasking
% 3. Create hemispheres
% 4. Save mutually exclusive masks
% 5. Create tsv-sidecar containing the names of the ROIs
% 6. Visual QC
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_im_Lesion2Mask('/MyStudy/sub-001/Lesion_T1_1.nii.gz', x);
% __________________________________
% Copyright 2015-2020 ExploreASL


%% --------------------------------------------------------
%% Admin
if nargin<1
    error('Missing input arguments');
elseif ~xASL_exist(LesionPath,'file')
    fprintf('%s\n',['Skipped because mask didnt exist: ' LesionPath]);
	return;
end

% Distinguish between lesion & ROI masks
[Fpath, Ffile] = xASL_fileparts(LesionPath);
if ~isempty(regexp(Ffile,'ROI')) && nargin>1 && isfield(x.D, 'ROICheckDir')
	fprintf('%s\n','Creating ROI masks');
    ImageSaveDir = x.D.ROICheckDir;
elseif ~isempty(regexp(Ffile,'Lesion')) && nargin>1 && isfield(x.D, 'LesionCheckDir')
	fprintf('%s\n','Creating Lesion masks');
	ImageSaveDir = x.D.LesionCheckDir;
else
	warning('Wrong lesion/ROI mask name');
	ImageSaveDir = Fpath;
end



%% --------------------------------------------------------
%% 1. If lesion is empty, skip this & delete the file
LesionIM = xASL_io_Nifti2Im(LesionPath);

if sum(LesionIM(:))<0.00001
    clear LesionIM
    xASL_delete(LesionPath);
    warning([LesionPath ' removed because it was empty']);
    return;
end

fprintf('%s\n','Printing lesion QC image');
LesionIM = xASL_im_ConvertMap2Mask(LesionIM(:,:,:,1));

DistanceMap = xASL_im_DistanceTransform(LesionIM);
DistanceMap = 1.5 .* DistanceMap; % convert voxels to 1.5 mm
PeriMask = DistanceMap>0 & DistanceMap<=16.6; % 25 mm
PeriMask(LesionIM) = 0; % make mutually exclusive


%% 2. BrainMasking
pGM         = xASL_io_Nifti2Im(x.P.Pop_Path_rc1T1);
pWM         = xASL_io_Nifti2Im(x.P.Pop_Path_rc2T1);
BrainMask   = (pGM+pWM)>0.5;
PeriMask    = PeriMask.*BrainMask;
AddMask     = LesionIM+PeriMask;
ContraMask  = xASL_im_Flip(LesionIM,1).*BrainMask;
ContraPeri  = xASL_im_Flip(PeriMask,1).*BrainMask;


%% 3. Create hemispheres
LeftHemisphere  = ones(size(BrainMask));
RightHemisphere = ones(size(BrainMask));
LeftHemisphere(1:round(size(LeftHemisphere,1)/2),:,:)     = 0;
RightHemisphere(round(size(RightHemisphere,1)/2):end,:,:) = 0;

if      sum(sum(sum(LeftHemisphere(logical(AddMask)))))>sum(sum(sum(RightHemisphere(logical(AddMask)))))
        Hemisphere      = (LeftHemisphere .* BrainMask) & ~AddMask;
        ContraLateral   = RightHemisphere .* BrainMask;
else    Hemisphere      = (RightHemisphere .* BrainMask) & ~AddMask;
        ContraLateral   = LeftHemisphere .* BrainMask;
end


%% 4. Save mutually exclusive masks
% 1) Intratumoral
LesionIM = uint8(LesionIM);
% 2) Peri, pGM+pWM
LesionIM(logical(PeriMask)) = 2;
% 3) Hemisphere
Hemisphere(logical(LesionIM)) = 0;
LesionIM(Hemisphere) = 3;
% 4) Contralateral intratumoral
LesionIM(logical(ContraMask)) = 4;
% 5) Contralateral peritumoral
LesionIM(logical(ContraPeri)) = 5;
% 6) Contralateral hemisphere
ContraLateral(logical(ContraMask) | logical(ContraPeri)) = 0;
LesionIM(logical(ContraLateral)) = 6;

xASL_io_SaveNifti(LesionPath, LesionPath, LesionIM, 8);


%% 5. Create tsv-sidecar containing the names of the ROIs
ROInames = {'Intratumoral' 'Perimask' 'Ipsilateral_Hemisphere' 'Contralateral_Intratumoral' 'Contralateral_Perimask' 'Contralateral_Hemisphere'};
PathTSV = fullfile(x.D.PopDir, [Ffile '.tsv']);
xASL_tsvWrite(ROInames, PathTSV, 1);


%% 6. Visual QC
% First create the individual mask overlays with different colors
for iIm=1:6
    OutIm{iIm} = xASL_vis_CreateVisualFig(x, {x.P.Pop_Path_rT1 LesionIM==iIm}, [], [0.75 0.35], [], {x.S.gray x.S.colors_ROI{iIm}});
end

% Initialize final image
OutImFinal = zeros(size(OutIm{1}));
OutImFinal(OutIm{1}==OutIm{2}) = OutImFinal(OutIm{1}==OutIm{2});

for iIm=1:length(OutIm)
    % Create masks with non-grey voxels (i.e. the overlay)
    OutImMask{iIm} = OutIm{iIm}(:,:,1)~=OutIm{iIm}(:,:,2) | OutIm{iIm}(:,:,2)~=OutIm{iIm}(:,:,3)  | OutIm{iIm}(:,:,1)~=OutIm{iIm}(:,:,3);
    OutImMask{iIm} = repmat(OutImMask{iIm}, [1 1 3]);
    
    % Add the overlay from each mask
    OutImFinal(OutImMask{iIm}) = OutIm{iIm}(OutImMask{iIm});
end

PathJPG = fullfile(ImageSaveDir, ['Regions_' Ffile '.jpg']);
xASL_adm_CreateDir(ImageSaveDir);
xASL_vis_Imwrite(OutImFinal, PathJPG);

    
end