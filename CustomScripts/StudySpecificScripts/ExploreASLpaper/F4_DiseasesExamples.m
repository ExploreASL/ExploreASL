x = ExploreASL_Master('',0);

%% Visualization admin
x.S.TraSlices = 53; % 35]; % Z slice
x.S.CorSlices = 75; %
x.S.SagSlices = 70; % left right
x.S.ConcatSliceDims = 1;
x.S.bCrop = 0;
% x.S.Square = 0;
bWhite = false;

clear ROOT IM
ROOT{1}     = 'C:\BackupWork\ASL\Novice\analysis\Population\Templates'; % diseased kids
ROOT{2}     = 'C:\BackupWork\ASL\EPAD_xASL_Paper\analysis_Sleep\Population\Templates'; % healthy adults
ROOT{3}     = 'C:\BackupWork\ASL\EPAD_xASL_Paper\analysis_EPAD\Population\Templates'; % diseased elderly

Modalities = {{'T1_bs-mean_Unmasked' 'pGM_bs-mean' 'pWM_bs-mean'} {'FLAIR_bs-mean_Unmasked' 'WMH_SEGM_bs-mean'}...
            {'CBF_bs-mean'} {'CBF_bs-CoV'} {'SD_bs-mean'} {'SNR_bs-mean'} {'mean_control_bs-mean'} {'noSmooth_M0_bs-mean'}...
            {'T1_bs-mean_Unmasked' 'SliceGradient_bs-mean_Unmasked'}};


ColorScale = {{x.S.gray x.S.red x.S.blue} {x.S.gray x.S.jet256} x.S.jet256 x.S.jet256 x.S.jet256...
            x.S.jet256 x.S.gray x.S.gray {x.S.gray x.S.jet256}};
Intens = {[0.75 0.2 0.2] [0.75 0.75] 1 1 1 1 1 1 [0.5 0.35]};
MaxClip = {[Inf Inf Inf] [Inf 0.7] 100 100 80 5 5000 7500000 [Inf 25]};

%% Create brainmask
pGM = xASL_io_Nifti2Im(fullfile(ROOT{2},'pGM_bs-mean.nii'));
pWM = xASL_io_Nifti2Im(fullfile(ROOT{2},'pWM_bs-mean.nii'));
pCSF = xASL_io_Nifti2Im(fullfile(ROOT{2},'pCSF_bs-mean.nii'));
bMask = (pGM+pWM+pCSF)>0.5;
bMasks3 = {bMask bMask bMask};
ZeroMasks3 = {ones(121,145,121) ones(121,145,121) ones(121,145,121)};
MasksAre = {ZeroMasks3 ZeroMasks3 bMasks3 bMasks3 bMasks3 bMasks3 ZeroMasks3 ZeroMasks3 ZeroMasks3};
EmptyWindow = {[] [] []};
MaxWindow = {EmptyWindow EmptyWindow EmptyWindow EmptyWindow EmptyWindow EmptyWindow...
            EmptyWindow EmptyWindow {[] 25}};

% 1) T1 + pGM+pWM
% 2) FLAIR +pWMH
% 3) CBF
% 4) bsCoV
% 5) tSD
% 6) tSNR
% 7) mean_control
% 8) M0
% 9) SliceGradient

%% Create the Figure
clear TotalIM TotalTotalIM
for iScanType=1:length(Modalities)
    clear IM Fpath
    for iStudy=1:length(ROOT)
        ProcessThis = true;
        clear tIM
        for iLayer = 1:length(Modalities{iScanType})
            Fpath{iScanType}{iLayer} = fullfile(ROOT{iStudy}, [Modalities{iScanType}{iLayer} '.nii']);
            if ~xASL_exist(Fpath{iScanType}{iLayer})
                ProcessThis = false;
            else % clipping
                tIM{iLayer} = xASL_io_Nifti2Im(Fpath{iScanType}{iLayer});
                tIM{iLayer}(tIM{iLayer}>MaxClip{iScanType}(iLayer)) = MaxClip{iScanType}(iLayer);
            end
        end
        if ProcessThis
            IM{iStudy} = xASL_im_CreateVisualFig(x, tIM, [], Intens{iScanType}, [], ColorScale{iScanType},false,MasksAre{iScanType},bWhite, MaxWindow{iScanType});
            % figure(1);imshow(IM,'InitialMagnification',250)
        elseif bWhite
            IM{iStudy} = ones(122,326,3);
        else
            IM{iStudy} = zeros(122,326,3);
        end
    end
    TotalIM{iScanType} = [IM{1} IM{2} IM{3}];
    figure(1);imshow(TotalIM{iScanType},'InitialMagnification',250)
end

TotalTotalIM = TotalIM{1};
for iTotal=2:length(TotalIM)
    TotalTotalIM = [TotalTotalIM;TotalIM{iTotal}];
end

figure(1);imshow(TotalTotalIM,'InitialMagnification',250) % export as EPS CMYK

%% Create the colorbars
for iBar=1:length(ColorScale)
    MaxValue = MaxClip{iBar}(end);
    if isfinite(MaxValue)
        DummyIm = repmat([0:0.01:1].*MaxValue,[101 1]);
        ColorMap = ColorScale{iBar};
        if length(ColorMap)<64
            ColorMap = ColorMap{end};
        end
        figure(iBar); imshow(DummyIm,[],'colormap',ColorMap,'InitialMagnification',400);
        colorbar;
    end
end
