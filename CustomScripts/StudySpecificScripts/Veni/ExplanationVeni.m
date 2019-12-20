%% Explanation ABBA

clear indIM

SDim            = xASL_io_Nifti2Im('C:\Backup\ASL\GENFI\GENFI_DF2_2_Long\analysis\dartel\Templates_Philips_Bsup\Template_sd_CBF.nii');
meanIM          = xASL_io_Nifti2Im('C:\Backup\ASL\GENFI\GENFI_DF2_2_Long\analysis\dartel\Templates_Philips_Bsup\Template_mean_CBF.nii');

indIMdir        = 'C:\Backup\ASL\SleepStudy\Analysis3\dartel';

Flist           = xASL_adm_GetFileList(indIMdir,'^qCBF_untreated.*\.(nii|nii\.gz)$');

for iM=1:18
    indIM(:,:,:,iM)     = xASL_io_Nifti2Im(Flist{iM});
end

indIM           = xASL_stat_MeanNan(indIM,4);

% dip_image(indIM)

IndExampleFile  = 'C:\Backup\ASL\Veni\Presentation\Normal_Example.nii';
Zmap_File       = 'C:\Backup\ASL\Veni\Presentation\Zmap_Healthy.nii';

pGM             = xASL_io_Nifti2Im(fullfile(x.D.MapsDir,'rgrey.nii'))>0.5;

% xASL_io_SaveNifti(x.D.ResliceRef,IndExampleFile,indIM);

% indIM           = indIM;
% meanIM          = meanIM;
% SDim            = SDim;

RatioIM         = (indIM./meanIM).*x.skull;
RatioMean       = median(RatioIM(pGM & isfinite(RatioIM)));
indIM           = indIM./RatioMean;

indIM           = xASL_im_ndnanfilter(indIM,'gauss',[1.885 1.885 1.885]);

Zmap            = xASL_im_ndnanfilter((x.skull.*(indIM-meanIM))./SDim,'gauss',[1.885 1.885 1.885]);

xASL_io_SaveNifti(x.D.ResliceRef,Zmap_File,Zmap);

Zmap(isnan(Zmap))       = 0;

jet_256                 = jet(256);
jet_256(129,:)          = 0;

Zmap                        = xASL_im_ndnanfilter(Zmap,'gauss',[1.885 1.885 1.885]);
Zmap(Zmap<0.1 & Zmap>-0.1)  = 0.5;
% Zmap                        = xASL_im_ndnanfilter(Zmap,'gauss',[1.885 1.885 1.885]);

figure(1);imshow(imrotate(Zmap(:,:,53).*x.skull(:,:,53),90),[-3 3],'colormap',jet_256,'border','tight','InitialMagnification',250)

dip_image([indIM meanIM])


%% Mean+Vasc from PreDiva

clear SDim meanIM
RatioFactor     = 3;

meanIM(:,:,:,1)   = xASL_io_Nifti2Im('C:\Backup\ASL\PreDiva\analysis\dartel\Templates\Template_mean_CBF_untreated_Ingenia.nii');
meanIM(:,:,:,2)   = xASL_io_Nifti2Im('C:\Backup\ASL\PreDiva\analysis\dartel\Templates\Template_mean_CBF_untreated_Intera.nii');
meanIM(:,:,:,3)   = RatioFactor.*xASL_io_Nifti2Im('C:\Backup\ASL\PreDiva\analysis\dartel\Templates\Template_mean_SD_Ingenia.nii');
meanIM(:,:,:,4)   = RatioFactor.*xASL_io_Nifti2Im('C:\Backup\ASL\PreDiva\analysis\dartel\Templates\Template_mean_SD_Intera.nii');
meanIM(:,:,:,5)   = xASL_io_Nifti2Im('C:\ExploreASL\Maps\Templates\MaxVesselTemplate.nii');

SDim(:,:,:,1)     = xASL_io_Nifti2Im('C:\Backup\ASL\PreDiva\analysis\dartel\Templates\Template_sd_CBF_untreated_Ingenia.nii');
SDim(:,:,:,2)     = xASL_io_Nifti2Im('C:\Backup\ASL\PreDiva\analysis\dartel\Templates\Template_sd_CBF_untreated_Intera.nii');
SDim(:,:,:,3)     = RatioFactor.*xASL_io_Nifti2Im('C:\Backup\ASL\PreDiva\analysis\dartel\Templates\Template_sd_SD_Ingenia.nii');
SDim(:,:,:,4)     = RatioFactor.*xASL_io_Nifti2Im('C:\Backup\ASL\PreDiva\analysis\dartel\Templates\Template_sd_SD_Intera.nii');

SDim            = xASL_stat_MeanNan(SDim,4);

SDim2           = (SDim.*(xASL_io_Nifti2Im('C:\ExploreASL\Maps\Templates\MRA_Viviani_xASL.nii').^0.5)).^0.7;

SDim            = SDim+SDim2.*9;

meanIM          = xASL_stat_MeanNan(meanIM,4);


% dip_image(meanIM)
% dip_image(SDim)

VascMeanFile    = 'C:\Backup\ASL\Veni\Presentation\VascMeanFile.nii';
VascSDFile      = 'C:\Backup\ASL\Veni\Presentation\VascSDFile.nii';

xASL_io_SaveNifti(x.D.ResliceRef,VascMeanFile,meanIM);
xASL_io_SaveNifti(x.D.ResliceRef,VascSDFile,SDim);

%% Explanation ABBA

clear indIM

SDim            = xASL_io_Nifti2Im('C:\Backup\ASL\GENFI\GENFI_DF2_2_Long\analysis\dartel\Templates_Philips_Bsup\Template_sd_CBF.nii');
meanIM          = xASL_io_Nifti2Im('C:\Backup\ASL\GENFI\GENFI_DF2_2_Long\analysis\dartel\Templates_Philips_Bsup\Template_mean_CBF.nii');
indIM           = xASL_io_Nifti2Im('C:\Backup\ASL\Veni\Presentation\VeniCBFExample2.nii');

% dip_image(indIM)

Zmap_File       = 'C:\Backup\ASL\Veni\Presentation\Zmap_Vasc.nii';

pGM             = xASL_io_Nifti2Im(fullfile(x.D.MapsDir,'rgrey.nii'))>0.5;

% xASL_io_SaveNifti(x.D.ResliceRef,IndExampleFile,indIM);

% indIM           = indIM;
% meanIM          = meanIM;
% SDim            = SDim;

RatioIM         = (indIM./meanIM).*x.skull;
RatioMean       = median(RatioIM(pGM & isfinite(RatioIM)));
indIM           = indIM./RatioMean;

indIM           = xASL_im_ndnanfilter(indIM,'gauss',[1.885 1.885 1.885]);

Zmap            = xASL_im_ndnanfilter((x.skull.*(indIM-meanIM))./SDim,'gauss',[1.885 1.885 1.885]);

xASL_io_SaveNifti(x.D.ResliceRef,Zmap_File,Zmap);

Zmap(isnan(Zmap))       = 0;

jet_256                 = jet(256);
jet_256(129,:)          = 0;

Zmap                        = xASL_im_ndnanfilter(Zmap,'gauss',[1.885 1.885 1.885]);
Zmap(Zmap<0.1 & Zmap>-0.1)  = 0.5;
% Zmap                        = xASL_im_ndnanfilter(Zmap,'gauss',[1.885 1.885 1.885]);

figure(1);imshow(imrotate(Zmap(:,:,53).*x.skull(:,:,53),90),[-3 3],'colormap',jet_256,'border','tight','InitialMagnification',250)

dip_image([indIM meanIM])

%% Explanation ABBA with correct Vasc

clear indIM SDim meanIM Zmap pGM RatioIM RatioMean

SDim            = xASL_io_Nifti2Im('C:\Backup\ASL\Veni\Presentation\VascSDFile.nii');
meanIM          = xASL_io_Nifti2Im('C:\Backup\ASL\Veni\Presentation\VascMeanFile.nii');
indIM           = xASL_io_Nifti2Im('C:\Backup\ASL\Veni\Presentation\VeniCBFExample2.nii');

% dip_image(indIM)

Zmap_File       = 'C:\Backup\ASL\Veni\Presentation\Zmap_Vasc2.nii';

pGM             = xASL_io_Nifti2Im(fullfile(x.D.MapsDir,'rgrey.nii'))>0.5;

% xASL_io_SaveNifti(x.D.ResliceRef,IndExampleFile,indIM);

% indIM           = indIM;
% meanIM          = meanIM;
% SDim            = SDim;

RatioIM         = (indIM./meanIM).*x.skull;
RatioMean       = median(RatioIM(pGM & isfinite(RatioIM)));
indIM           = indIM./RatioMean;

indIM           = xASL_im_ndnanfilter(indIM,'gauss',[1.885 1.885 1.885]);

Zmap            = xASL_im_ndnanfilter((x.skull.*(indIM-meanIM))./SDim,'gauss',[1.885 1.885 1.885]);

xASL_io_SaveNifti(x.D.ResliceRef,Zmap_File,Zmap);

Zmap(isnan(Zmap))       = 0;

jet_256                 = jet(256);
jet_256(129,:)          = 0;

Zmap                        = xASL_im_ndnanfilter(Zmap,'gauss',[1.885 1.885 1.885]);
Zmap(Zmap<0.1 & Zmap>-0.1)  = 0.5;
% Zmap                        = xASL_im_ndnanfilter(Zmap,'gauss',[1.885 1.885 1.885]);

figure(1);imshow(imrotate(Zmap(:,:,53).*x.skull(:,:,53),90),[-3 3],'colormap',jet_256,'border','tight','InitialMagnification',250)

dip_image([indIM meanIM])
