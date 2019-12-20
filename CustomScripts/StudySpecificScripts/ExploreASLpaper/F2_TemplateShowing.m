%% This shows the different templates that are available through ExploreASL

ExploreASL_Master('',0);
clear
%% TemplateShowing for ExploreASL poster
%  for the paper this was done with mricron stacking the transversal slices side to side

TemplatesDir    = 'C:\Users\kyrav\Desktop\Gdrive\XploreLab\ProjectsPending\ExploreASL manuscript\Figures\Fig2_Templates\SequenceAverages';

pGM             = xASL_io_Nifti2Im('C:\ExploreASL\Maps\rgrey.nii');
pGM             = pGM>0.5;
bMASK           = xASL_io_Nifti2Im('C:\ExploreASL\Maps\rbrainmask.nii');

%%

PathName{1}        = fullfile(TemplatesDir, '1_GE_3Dspiral_WIP_CBF_INOT.nii'); % Fernando INOX (alternative Joost Kuijer, but takes lots of time to get AD M0 data...)
PathName{2}        = fullfile(TemplatesDir, '2_GE_3D_spiral_Product_CBF_VESPA.nii'); % VESPA (previously GENFI 3D spiral (from registration paper)
PathName{3}        = fullfile(TemplatesDir, '3_Philips_2D_EPI_PCASL_CBF_Sleep.nii'); % Sleep data Oslo 2015 (previously GENFI data)
PathName{4}        = fullfile(TemplatesDir, '4_Philips_3DGRASE_1800ms_PLD_CBF_Philips.nii'); % Kim vd Ven, Philips
PathName{5}        = fullfile(TemplatesDir, '5_Siemens_2D_EPI_PCASL_noBsup_CBF_AAS.nii'); % BoleStudien (alternative Chris Chen, alternative BioCog) -> need to rerun for without masking
PathName{6}        = fullfile(TemplatesDir, '6_Siemens_3DGRASE_PCASL_CBF_MoodStudy.nii'); % Mood (was GENFI (alternative Mood study Rik)

%%
clear IM S

jet_256         = jet(256);
jet_256(1,:)    = 0;

for iN=1:length(PathName)
    tIM             = xASL_io_Nifti2Im(PathName{iN});
    meanIM(iN)      = mean(tIM(isfinite(tIM) & pGM));
    IM(iN,:,:,:)    = tIM.*60./meanIM(iN).*bMASK;
end

x.S.TraSlices     = [53 73];
x.S.CorSlices     = [75 98];
x.S.SagSlices     = [64 86];
x.S.ConcatSliceDims = 1;
x.S.Square        = 0;

DATA_OUT    = xASL_im_TransformData2View( IM, x );
TotalIM     = xASL_im_TileImages(DATA_OUT,1);
figure(1);imshow(TotalIM,[0 125]) % ,'colormap',jet_256

%% --------------------------------------------------------------------
%% Source images showing for ExploreASL paper

clear Fname

% do as if healthy controls, as you don't really see the difference between healthy & diseased, as much as you see age differences
% in these average sequence images.
% Specifically: vascular peaks are smoothed away in 3D readout images
% Specifically: atrophy is far reduced by warping to MNI
% CBF intensity differences are removed by scaling



PathName{1}       = fullfile(TemplatesDir, '1_GE_3Dspiral_WIP_M0_INOT.nii');
PathName{2}       = fullfile(TemplatesDir, '2_GE_3D_spiral_Product_M0_VESPA.nii');
PathName{3}       = fullfile(TemplatesDir, '3_Philips_2D_EPI_PCASL_Control_Sleep.nii');
PathName{4}       = fullfile(TemplatesDir, '4_Philips_3DGRASE_1800ms_PLD_Control_Philips.nii');
PathName{5}       = fullfile(TemplatesDir, '5_Siemens_2D_EPI_PCASL_noBsup_Control.nii');
PathName{6}       = fullfile(TemplatesDir, '6_Siemens_3DGRASE_PCASL_Control_MoodStudy.nii');


%%
clear IM x.S meanIM

jet_256         = jet(256);
jet_256(1,:)    = 0;

for iN=1:length(PathName)
    tIM                 = xASL_io_Nifti2Im(PathName{iN});
    meanIM(iN)          = mean(tIM(isfinite(tIM) & pGM));
    IMmask(iN,:,:,:)    = tIM.*60./meanIM(iN).*bMASK;
    IM(iN,:,:,:)        = tIM.*60./meanIM(iN);
end

% Correction window/level:
IM(3,:,:,:)            = IM(3,:,:,:).*0.75;

x.S.TraSlices     = 53; % [53 73];
x.S.CorSlices     = 75; % [75 98];
x.S.SagSlices     = 64; % [64 86];
x.S.ConcatSliceDims = 1;
x.S.Square        = 0;

DATA_OUTmask    = xASL_im_TransformData2View( IMmask, x );
DATA_OUT        = xASL_im_TransformData2View( IM, x );

TotalIMmask     = xASL_im_TileImages(DATA_OUTmask,1);
TotalIM         = xASL_im_TileImages(DATA_OUT,1);
figure(2);imshow(TotalIMmask,[0 125],'colormap',grey,'border','tight') % masked
figure(3);imshow(TotalIM,[0 125],'colormap',grey,'border','tight') % without masking
