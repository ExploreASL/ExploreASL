%% Create images for movie
% [x]     = InitializeExploreASL;

clear PseudoCBFim

SmoothBrainMask     = xASL_io_Nifti2Im('C:\ExploreASL\Maps\rbrainmask.nii');

Vasc_IM     = xASL_io_Nifti2Im('C:\ExploreASL\Maps\Templates\MaxVesselTemplate.nii');
Bias_IM     = xASL_io_Nifti2Im('C:\ExploreASL\Maps\Templates\ATT_BiasField.nii');
Mean_IM     = xASL_io_Nifti2Im('C:\ExploreASL\Maps\Templates\Philips_2DEPI_Bsup_CBF_Template.nii');

% x.S.TraSlices     = [53 62 74 87];
% x.S.CorSlices     = [53 62 74 87];
% x.S.SagSlices     = [53 62 74 87];
% x.S.ConcatSliceDims = 0;

% Vasc_IM             = TransformDataViewDimension(Vasc_IM,[],S);
% Bias_IM             = TransformDataViewDimension(Bias_IM,[],S);
% Mean_IM             = TransformDataViewDimension(Mean_IM,[],S);
% SmoothBrainMask     = TransformDataViewDimension(SmoothBrainMask,[],S);

Vasc_IM             = imrotate(Vasc_IM(:,:,53),90);
Bias_IM             = imrotate(Bias_IM(:,:,53),90);
Mean_IM             = imrotate(Mean_IM(:,:,53),90);
SmoothBrainMask     = imrotate(SmoothBrainMask(:,:,53),90);

spatCoV     = fliplr([0.1:0.05:0.45 0.5:0.025:2.5]);
for iS=1:length(spatCoV)
    PseudoCBFim(:,:,iS) = (Mean_IM - Bias_IM.*5.*spatCoV(iS).^2 + Vasc_IM./5.*spatCoV(iS).^2) .* SmoothBrainMask;
end

SmoothBrainMask     = logical(repmat(SmoothBrainMask,[1 1 size(PseudoCBFim,3)]));

PseudoCBFim(isnan(PseudoCBFim))     = 0;
PseudoCBFim(SmoothBrainMask)        = max(0,PseudoCBFim(SmoothBrainMask));
% PseudoCBFim(PseudoCBFim==0)         = 105;


%% Save movie as 3D NIfTI & upscale to visible movie resolution
RootDir     = 'C:\Users\henkj\Dropbox\Itinerant Science\Veni\Veni Interview\Images';
SaveFile    = fullfile(RootDir, 'PseudoCBFim.gif');
SaveNII     = fullfile(RootDir, 'PseudoCBFim.nii');

if  exist(SaveFile,'file'); delete(SaveFile);end
if  exist(SaveNII,'file');  delete(SaveNII) ;end

xASL_io_SaveNifti(x.D.ResliceRef, SaveNII, PseudoCBFim);
UpsampleNII( SaveNII, SaveNII, [0.5 0.5 1]);
PseudoCBFim     = xASL_io_Nifti2Im(SaveNII);

%% Adapt range to colorscale
jet_256     = double(jet(256));
jet_256(1,:)= 0;

MaxInt          = 100;
PseudoCBFim(PseudoCBFim<1)   = 1;
PseudoCBFim     = double(round(PseudoCBFim./MaxInt.*255));
PseudoCBFim(PseudoCBFim>255)    = 255;
DiffMin         = min(PseudoCBFim(:))-1;
PseudoCBFim     = PseudoCBFim - DiffMin;

PseudoCBFim     = PseudoCBFim(:,:,[61:2:end]);

for iP=1:size(PseudoCBFim,3)
    ImIndex     =   PseudoCBFim(:,:,iP);
    if iP==1
        imwrite(ImIndex,jet_256,SaveFile,'gif', 'Loopcount',inf, 'DelayTime',0); %,'Screensize',[80 80]
    else
        imwrite(ImIndex,jet_256,SaveFile,'gif', 'WriteMode','append', 'DelayTime',0);
    end
end
