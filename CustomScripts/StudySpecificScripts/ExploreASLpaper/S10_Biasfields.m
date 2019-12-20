%% S10 biasfields

clear
x = ExploreASL_Master('',0);

%% TemplateShowing for ExploreASL poster
%  for the paper this was done with mricron stacking the transversal slices side to side

x.WBmaskNarrow = xASL_io_Nifti2Im(fullfile(x.D.MapsDir, 'WBmaskASLnarrow.nii'));

%% ------------------------------------------------------------------------------------------------------------
%% Admin

x.D.TemplatesStudyDir = 'C:\Users\kyrav\Desktop\Gdrive\XploreLab\ProjectsPending\ExploreASL manuscript\Figures\Fig2_Templates\SequenceAverages';

SiteAveragePath{1}        = fullfile(x.D.TemplatesStudyDir, '1_GE_3Dspiral_WIP_CBF_INOT.nii'); % Fernando INOX (alternative Joost Kuijer, but takes lots of time to get AD M0 data...)
SiteAveragePath{2}        = fullfile(x.D.TemplatesStudyDir, '2_GE_3D_spiral_Product_CBF_VESPA.nii'); % VESPA (previously GENFI 3D spiral (from registration paper)
SiteAveragePath{3}        = fullfile(x.D.TemplatesStudyDir, '3_Philips_2D_EPI_PCASL_CBF_Sleep.nii'); % Sleep data Oslo 2015 (previously GENFI data)
SiteAveragePath{4}        = fullfile(x.D.TemplatesStudyDir, '4_Philips_3DGRASE_1800ms_PLD_CBF_Philips.nii'); % Kim vd Ven, Philips
SiteAveragePath{5}        = fullfile(x.D.TemplatesStudyDir, '5_Siemens_2D_EPI_PCASL_noBsup_CBF_AAS.nii'); % BoleStudien (alternative Chris Chen, alternative BioCog) -> need to rerun for without masking
SiteAveragePath{6}        = fullfile(x.D.TemplatesStudyDir, '6_Siemens_3DGRASE_PCASL_CBF_MoodStudy.nii'); % Mood (was GENFI (alternative Mood study Rik)

%% Hack xASL_wrp_CreateBiasfield

pGMname = 'C:\ExploreASL\Maps\rgrey.nii';

SiteSet = 1;
x.S.SetsID(:,SiteSet) = [1:6];
AllSites = unique( x.S.SetsID(:,SiteSet) );
x.S.SetsOptions{SiteSet} = {'GEspiral' 'GEspiralWIP' 'Philips_2DEPI' 'Philips_3DGRASE' 'Siemens_2DEPI' 'Siemens_3DGRASE'};
SiteNames = x.S.SetsOptions{SiteSet};
nSites = length(AllSites);

pGM = xASL_io_Nifti2Im(pGMname)>0.5;


%% Import average site-CBF images
for iSite=1:nSites
%     SiteAveragePath{iSite} = fullfile(x.D.TemplatesStudyDir, [CBF_prefix '_Site_' SiteNames{iSite} '_bs-mean.nii']);

    if xASL_exist(SiteAveragePath{iSite},'file')
        RobustMean(:,:,:,iSite) = xASL_io_Nifti2Im(SiteAveragePath{iSite});
    end
end


%% COPY PASTE FROM xASL_wrp_CreateBiasfield



    % This needs improvement, to avoid PV effects in the intensity normalization
    % ------------------------------------------------------------------------------------------------------------
    % Create smooth images
    % Store reference images
    RefIM{1} = xASL_im_rotate(RobustMean(:,:,53,1),90);
    for iSite=2:nSites
        % 1) visualize mean CBF images
        RefIM{1} = [RefIM{1} xASL_im_rotate(RobustMean(:,:,53,iSite),90)];
    end
    
    MaskIM = repmat(x.WBmaskNarrow,[1 1 1 size(RobustMean,4)]);
    RobustMean(~MaskIM) = NaN;
    Ind_Biasfield = RobustMean;

    RefIM{2} = xASL_im_rotate(Ind_Biasfield(:,:,53,1),90);
    for iSite=2:nSites
        % 2) visualize mean CBF images after masking
        RefIM{2} = [RefIM{2} xASL_im_rotate(Ind_Biasfield(:,:,53,iSite),90)];
    end    
    
    fprintf('%s','Smoothing images:   ');
    for iSite=1:nSites
        xASL_TrackProgress(iSite,nSites);
        % Wide smoothing & expand
        Ind_Biasfield(:,:,:,iSite) = xASL_im_ndnanfilter(Ind_Biasfield(:,:,:,iSite),'gauss',[16 16 16],0);
        while sum(sum(sum(isnan(Ind_Biasfield(:,:,:,iSite)))))>0
              Ind_Biasfield(:,:,:,iSite) = xASL_im_ndnanfilter(Ind_Biasfield(:,:,:,iSite),'gauss',[8 8 8],2);
        end
    end
    fprintf('\n');

    RefIM{3} = xASL_im_rotate(Ind_Biasfield(:,:,53,1),90);
    for iSite=2:nSites
        % 3) visualize mean CBF images after masking & smoothing with extrapolation
        RefIM{3} = [RefIM{3} xASL_im_rotate(Ind_Biasfield(:,:,53,iSite),90)];
    end


    % ------------------------------------------------------------------------------------------------------------
    % Create average bias-field
    AV_biasfield = xASL_stat_MeanNan(Ind_Biasfield,4);
    % Turn NaNs into zeros
    AV_biasfield(isnan(AV_biasfield))= 0;
    % Turn biasfield to mean CBF of 60
    ScaleFactor = 60/xASL_stat_MeanNan(AV_biasfield(pGM));
    AV_biasfield = ScaleFactor .* AV_biasfield;

    % Put average images in RefIM
    % 4) average mean CBF
    RefIM{4} = repmat(xASL_im_rotate(AV_biasfield(:,:,53),90),[1 nSites]);





    % ------------------------------------------------------------------------------------------------------------
    % Create new bias-fields, masked & expand
    % using AV_biasfield = Ind_Biasfield .* BiasFieldMultipl + BiasFieldAdditiv
    BiasFieldMultipl = repmat(AV_biasfield,[1 1 1 nSites]) ./ Ind_Biasfield ; % .* repmat(Bmask,[1 1 1 nSites])
    BiasFieldMultipl(BiasFieldMultipl<=0) = NaN;

    fprintf('%s','Extrapolating biasfields:   ');
    for iSite=1:nSites
        xASL_TrackProgress(iSite,nSites);
        % Extrapolate values far outside brainmask to avoid division artifacts
        BiasFieldMultipl(:,:,:,iSite) = xASL_im_ndnanfilter( BiasFieldMultipl(:,:,:,iSite) ,'gauss',[14 14 14],2);
        BiasFieldMultipl(:,:,:,iSite) = xASL_im_ndnanfilter( BiasFieldMultipl(:,:,:,iSite) ,'gauss',[14 14 14],2);
    end
    
    fprintf('\n%s','Robust clipping:   ');
    for iSite=1:nSites
        xASL_TrackProgress(iSite,nSites);
        tempIM  = BiasFieldMultipl(:,:,:,iSite);
        for ii=1:8
            tempIM = xASL_im_ndnanfilter( tempIM ,'gauss',[1.885 1.885 1.885],0);
        end
        MaxThr = max(max(max(tempIM(pGM))));
        tempIM = BiasFieldMultipl(:,:,:,iSite);
        tempIM(tempIM>MaxThr) = MaxThr;
        BiasFieldMultipl(:,:,:,iSite) = tempIM;
    end




    % ------------------------------------------------------------------------------------------------------------
    % Store reference images
    RefIM{5} = xASL_im_rotate(BiasFieldMultipl(:,:,53,1),90);
    for iSite=2:nSites
        % 5) Biasfields
        RefIM{5} = [RefIM{5} xASL_im_rotate(BiasFieldMultipl(:,:,53,iSite),90)];
    end


%% Remove mean scaling between biasfields for visualization

RefVis = squeeze(xASL_im_rotate(BiasFieldMultipl(:,:,53,:),90));
ScaleF = squeeze(mean(mean(RefVis,1),2))';
for iSite=1:nSites
    RefVis(:,:,iSite) = RefVis(:,:,iSite)./ScaleF(iSite);
end

RefVis = [RefVis(:,:,1) RefVis(:,:,2) RefVis(:,:,3) RefVis(:,:,4) RefVis(:,:,5) RefVis(:,:,6)];
    
% TotalIM = [RefIM{1}; RefIM{2}; RefIM{3}; RefIM{4};; RefIM{5}.*30; RefIM{1}.*RefIM{5}];
TotalIM = [RefIM{1}; RefVis.*80; RefIM{1}.*RefIM{5}];
figure(1);imshow(TotalIM,[0 120],'colormap',x.S.jet256)
