%% Create ghosting template ExploreQC
% Here we want to distinguish the brainmask,
% and AP (and LR?) effects
% So we use diagonal lines to separate the AP (& LR?) areas outside the
% brainmask, and clip this with the most lateral brainlines
%
% These are the template ROIs for the Ghost-To-Signal Ratio
% REF: https://mriqc.readthedocs.io/en/stable/iqms/bold.html
% 1 = Signal
% 2 = Non-ghost
% 3 = Ghost

x = ExploreASL_Master('',0);

PathOutput = 'C:\ExploreASL\External\SPMmodified\MapsAdded\GhostSignalRatio.nii';

BrainmaskPath = 'C:\ExploreASL\External\SPMmodified\MapsAdded\brainmask.nii';
BrainmaskIm = xASL_io_Nifti2Im(BrainmaskPath);

% We first create a matrix where voxelvalues have distance from X-axis midline 
LRdistanceMatrix = zeros(size(BrainmaskIm));
APdistanceMatrix = zeros(size(BrainmaskIm));

DimIm = size(BrainmaskIm);
for iJ=1:DimIm(1)
    LRdistanceMatrix(iJ,:,:) = abs(iJ-(0.5*DimIm(1)));
end

% Then create a matrix where voxelvalues have distance from Y-axis midline
for iK=1:DimIm(2)
    APdistanceMatrix(:,iK,:) = abs(iK-(0.5*DimIm(2)));
end

% Now all voxels with same values are on the diagonals
% Voxels with APdistance>LRdistance are AP region
% Voxels with LRdistance>APdistance are lateral region

GhostImage = single(APdistanceMatrix>LRdistanceMatrix)+2;
GhostImage(logical(BrainmaskIm)) = 1;

% Everywhere it is 3 (anterior), 
% if we flip the brainmask,
% where brainmask=1, make GhostImage=2,
% to create the curves

% 1) take anterior part of brainmask
FlipMask = BrainmaskIm;
FlipMask(:,1:round(0.5*DimIm(2)),:) = 0;
FlipMask = FlipMask(:,[end:-1:1],:);
FlipMask(:,ceil(0.5*DimIm(2))+1:end,:) = FlipMask(:,1:floor(0.5*DimIm(2)),:);
FlipMask(:,1:floor(0.5*DimIm(2)),:) = 0;

GhostImageAnterior = GhostImage;
GhostImageAnterior(GhostImage==3 & FlipMask==0) = 2;

% 2% do the same for posterior part
FlipMask = BrainmaskIm;
FlipMask(:,round(0.5*size(FlipMask,2))+1:end,:) = 0;
FlipMask = FlipMask(:,[end:-1:1],:);
FlipMask(:,1:floor(0.5*DimIm(2)),:) = FlipMask(:,ceil(0.5*DimIm(2))+1:end,:);
FlipMask(:,floor(0.5*DimIm(2))+1:end,:) = 0;

GhostImagePosterior = GhostImage;
GhostImagePosterior(GhostImage==3 & FlipMask==0) = 2;

GhostImage(:,1:floor(0.5*DimIm(2)),:) = GhostImagePosterior(:,1:floor(0.5*DimIm(2)),:);
GhostImage(:,floor(0.5*DimIm(2))+1:end,:) = GhostImageAnterior(:,floor(0.5*DimIm(2))+1:end,:);

GhostImage(:,:,1:45) = 0; % remove inferior slices
GhostImage(:,:,86:end) = 0; % remove superior slices

xASL_io_SaveNifti(BrainmaskPath, PathOutput, uint8(GhostImage), [], 1);