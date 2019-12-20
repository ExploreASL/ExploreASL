%% Create volumetric ASL movie
% [x]     = InitializeExploreASL;

%% Create synthetic CBF image

clear PseudoCBFim

pGM     = xASL_io_Nifti2Im('C:\ExploreASL\Maps\Templates\rc1T1_ASL_res.nii');
pWM     = xASL_io_Nifti2Im('C:\ExploreASL\Maps\Templates\rc2T1_ASL_res.nii');
bMask   = xASL_io_Nifti2Im('C:\ExploreASL\Maps\Templates\brainmask.nii');
% CBF     = xASL_io_Nifti2Im('C:\ExploreASL\Maps\Templates\Philips_2DEPI_Bsup_CBF_Template.nii');

pGM     = pGM(:,:,53);
pWM     = pWM(:,:,53);
pWM(:)  = 0;
bMask   = logical(bMask(:,:,53));

pCSF    = bMask-pGM-pWM;

IM1     = pGM.*60+pWM.*20;
IM1     = IM1./max(IM1(:));
IM1     = xASL_im_rotate(IM1,90);

% pGM_small   = imerode(pGM,strel('disk',1));

SpatVol = [1:0.02:4];

for iS=1:length(spatCoV)
    pGM_frame           = pGM.^SpatVol(iS);
    pWM_frame           = max(0,pWM-pCSF.^(1/SpatVol(iS).^12));
    PseudoCBFim(:,:,iS) = max(0,pGM_frame+pWM_frame./3);
    PseudoCBFim(:,:,iS) = PseudoCBFim(:,:,iS) .* bMask;
end

PseudoCBFim             = xASL_im_rotate(PseudoCBFim,90);

%% Save movie as 3D NIfTI & upscale to visible movie resolution
RootDir     = 'C:\Users\henkj\Dropbox\Itinerant Science\Veni\Veni Interview\Images';
SaveFile    = fullfile(RootDir, 'PseudoCBF_struct_im.gif');
SaveNII     = fullfile(RootDir, 'PseudoCBF_struct_im.nii');

if  exist(SaveFile,'file'); delete(SaveFile);end
if  exist(SaveNII,'file');  delete(SaveNII) ;end

PseudoCBFim     = PseudoCBFim.^0.67;

xASL_io_SaveNifti(x.D.ResliceRef, SaveNII, PseudoCBFim);
UpsampleNII( SaveNII, SaveNII, [0.5 0.5 1]);
PseudoCBFim     = xASL_io_Nifti2Im(SaveNII);

%% Adapt range to colorscale
jet_256     = double(jet(256));
jet_256(1,:)= 0;

MaxInt          = 1;
PseudoCBFim(PseudoCBFim<0)   = 0;
PseudoCBFim     = double(round(PseudoCBFim./MaxInt.*255));
PseudoCBFim(PseudoCBFim>255)    = 255;
DiffMin         = min(PseudoCBFim(:))-1;
PseudoCBFim     = PseudoCBFim - DiffMin;

% PseudoCBFim     = PseudoCBFim(:,:,[2:3:end]);
PseudoCBFim     = PseudoCBFim(:,:,[40:3:end]);

for iP=1:size(PseudoCBFim,3)
    ImIndex     =   PseudoCBFim(:,:,iP);
    if iP==1
        imwrite(ImIndex,jet_256,SaveFile,'gif', 'Loopcount',inf, 'DelayTime',0); %,'Screensize',[80 80]
    else
        imwrite(ImIndex,jet_256,SaveFile,'gif', 'WriteMode','append', 'DelayTime',0);
    end
end
