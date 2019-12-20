%% Veni example/pilot data

x.D.ROOT    = 'C:\Backup\ASL\BioCog\DataFreeze1';
dataparfile     = fullfile(x.D.ROOT,'DATA_PAR_Berlin_exclude.m');
x         = InitializeExploreASL( dataparfile );

%% Load data
ModalityN       = {'qCBF_untreated' 'PV_pGM' 'PV_pWM'};
SuffixN         = {'_ASL_1' '' ''};
for iM=1 %1:length(ModalityN)
    % pre-allocate
    IMim{iM}    = zeros(x.nSubjects,145,121,121);
    for iS=1:x.nSubjects
        xASL_TrackProgress(iS,x.nSubjects);

        Fname                   = fullfile(x.D.PopDir,[ModalityN{iM} '_' x.SUBJECTS{iS} SuffixN{iM} '.nii']);
        if     ~exist(Fname,'file')
                Fname           = fullfile(x.D.PopDir,[ModalityN{iM} '_' x.SUBJECTS{iS} SuffixN{iM} '.nii.gz']);
        end

         tIM    = xASL_io_Nifti2Im( Fname );
%          IMim{iM}(iS,:,:,:)     = xASL_im_rotate(dip_array(gaussf(tIM,3)),90);
         IMim{iM}(iS,:,:,:)     = xASL_im_rotate(tIM,90);
    end
end

for ii=1 %1:3
    xASL_TrackProgress(ii,3);
    meanIM{ii}      = squeeze(xASL_stat_MeanNan(IMim{ii},1));
    stdIM{ii}       = squeeze(xASL_stat_StdNan(IMim{ii},1));
end

bMask       = xASL_im_rotate(x.SteepSkullMap,90);
% GMmask      = meanIM{2}>0.5;
% WMmask      = meanIM{3}>0.5;


% PV_GM       = meanIM{1}./meanIM{2};
% PV_WM       = meanIM{1}./meanIM{3};
% GM_CBF      = xASL_stat_SumNan(xASL_stat_SumNan(xASL_stat_SumNan(PV_GM(GMmask))))./sum(GMmask(:))
% WM_CBF      = xASL_stat_SumNan(xASL_stat_SumNan(xASL_stat_SumNan(PV_WM(WMmask))))./sum(WMmask(:))

% perform PVEC on group average
iP  =1;
for ii=1:10
    IND                 = squeeze(IMim{iP}(ii,:,:,53)).*bMask(:,:,53);
    Z_im                = (IND-meanIM{iP}(:,:,53)) .*bMask(:,:,53) ./stdIM{iP}(:,:,53).*bMask(:,:,53);
    TotalIM(:,:,ii)     = [meanIM{iP}(:,:,53).*bMask(:,:,53);IND.*bMask(:,:,53);abs(IND-meanIM{iP}(:,:,53)).*bMask(:,:,53);Z_im.*bMask(:,:,53)./10];
end

jet_256         = jet(256);
jet_256(1,:)    = 0;
TotalIM(TotalIM==0)     = NaN;
figure(1);imshow(singlesequencesort(TotalIM,10),[0 100],'Colormap',jet_256)
