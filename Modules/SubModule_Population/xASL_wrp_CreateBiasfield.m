function xASL_wrp_CreateBiasfield(x)
%xASL_wrp_CreateBiasfield Create sequence-specific intensity biasfields
%
% FORMAT: xASL_wrp_CreateBiasfield(x)
%
% INPUT:
%   x           - struct containing pipeline environment parameters
% OUTPUT: n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function creates a smooth biasfield as intensity map to
%              normalize intensities between different
%              sequences/scanner/sites within a single study.
%              This is a simple pragmatic approach and is not validated,
%              but is the best we can do.
%
%              First acquires average additive & multiplicative factors for total GM,
%              then does a smooth voxel-wise rescaling. This doesn't make an assumption
%              whether site or sequence differences are additive or multiplicative,
%              but rather fits them both. Global scaling it performed to GM CBF == 60 mL/100g/min
%
%              NB: make sure that sequence resolution differences have been
%              taken in account before creating these biasfields
%              PM: add normalization of between-subjects SD as well.
%              PM: are there other things we can normalize?
%
%              This function has the following sections:
%
%  0. Admin
%  1. Load images & store site-average CBF images
%  2. Smooth the site-average CBF images to biasfields
%  3. Create average biasfield
%  4. Create new bias-fields, masked & expand
%  5. Save biasfields
%  6. Create backup of each CBF image before applying biasfields
%  7. Print QC reference images
%  8. Rescale CBF images (apply the biasfield correction)
%  9. Re-load & store site-average CBF images
% 10. Perform site-wise SD correction & save new CBF images
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_wrp_CreateBiasfield(x);
% __________________________________
% Copyright 2015-2023 ExploreASL

    x.D.TemplatesStudyDir = fullfile(x.D.PopDir,'Templates');
    CBF_prefix = ['q' x.P.CBF];

    % Select the set that contains the variable Site.mat
    if isfield(x.S,'SetsName')
        iSiteSet = find(strcmp(x.S.SetsName,'Site'));
    end

    x.WBmaskNarrow = xASL_io_Nifti2Im(fullfile(x.D.MapsSPMmodifiedDir,'ParenchymNarrow.nii'));

    %% ------------------------------------------------------------------------------------------------------------
    %% 0. Admin
    pGMname = fullfile(x.D.MapsSPMmodifiedDir, 'rc1T1.nii');

    if ~exist('iSiteSet','var')
        fprintf('%s\n','No multiple sites found, BiasField creation and normalization skipped');
        return;
    end

    % Define sites
    AllSites = unique(x.S.SetsID(:,iSiteSet));
    SiteNames = x.S.SetsOptions{iSiteSet};
    nSites = numel(AllSites);

	% only run this function if there are multiple sites defined, & hasn't been run before    
    if ~(nSites>1 && isempty(xASL_adm_GetFileList(x.D.PopDir,'(?i)^Biasfield_.*_Site_\d*.*\.nii$')))
        fprintf('%s\n','No multiple sites found, BiasField creation and normalization skipped');
        return;
    end

    pGM = xASL_io_Nifti2Im(pGMname)>0.5;


    %% ------------------------------------------------------------------------------------------------------------
    %% Create average images per scanner/sequence/site
    % This should be specified in a variable "Site" column in
    % participants.tsv
    % If you have multiple sessions, then specify as how you want to
    % average/rescale them
    % NB: make sure you exclude the bad apples first, since these can mess up
    % this rescaling process (although robust mean average maps are created
    % here)

    fprintf('%s\n',['Creating biasfields for ' num2str(nSites) ' sites']);

    %% 1. Load images & store site-average CBF images
    for iSite=1:nSites
        IM = single.empty(0);
        SiteScans = find(x.S.SetsID(:,iSiteSet)==AllSites(iSite));
        nScans = numel(SiteScans);

        ScanN = 1;
        fprintf('%s', ['Loading scans for site ' num2str(iSite) ':    ']);
        for iScan=1:nScans
            xASL_TrackProgress(iScan,nScans);
            iSubSess = SiteScans(iScan);
            iSub = ceil(iSubSess/x.dataset.nSessions);
            iSess = iSubSess- ((iSub-1)*x.dataset.nSessions);
            FileName = fullfile(x.D.PopDir,[CBF_prefix '_' x.SUBJECTS{iSub} '_' x.SESSIONS{iSess} '.nii']);
            if xASL_exist(FileName, 'file')
                IM(:,ScanN) = xASL_im_IM2Column(xASL_io_Nifti2Im(FileName),x.S.masks.WBmask);
                ScanN = ScanN+1;
                NameList{ScanN,iSite} = [x.SUBJECTS{iSub} '_' x.SESSIONS{iSess}];
            end
        end

        NoOutliers = xASL_stat_RobustMean(IM); % use SoS function
        RobustMean(:,:,:,iSite) = xASL_im_Column2IM(xASL_stat_MeanNan(IM(:,NoOutliers), 2), x.S.masks.WBmask);
    end




    %% This needs improvement, to avoid PV effects in the intensity normalization
    %% ------------------------------------------------------------------------------------------------------------
    %% 2. Smooth the site-average CBF images to biasfields
    % Store reference images
    RefIM{1} = xASL_im_rotate(RobustMean(:,:,53,1),90);
    for iSite=2:nSites
        % Reference image 1: visualize mean CBF images
        RefIM{1} = [RefIM{1} xASL_im_rotate(RobustMean(:,:,53,iSite),90)];
    end

    MaskIM = repmat(x.WBmaskNarrow,[1 1 1 size(RobustMean,4)]);
    RobustMean(~MaskIM) = NaN;
    Ind_Biasfield = RobustMean;

    RefIM{2} = xASL_im_rotate(Ind_Biasfield(:,:,53,1),90);
    for iSite=2:nSites
        % Reference image 2: visualize mean CBF images after masking
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
        % Reference image 3: visualize mean CBF images after masking & smoothing with extrapolation
        RefIM{3} = [RefIM{3} xASL_im_rotate(Ind_Biasfield(:,:,53,iSite),90)];
    end


    %% ------------------------------------------------------------------------------------------------------------
    %% 3. Create average biasfield
    AV_biasfield = xASL_stat_MeanNan(Ind_Biasfield,4);
    % Turn NaNs into zeros
    AV_biasfield(isnan(AV_biasfield))= 0;
    % Turn biasfield to mean CBF of 60
    ScaleFactor = 60/xASL_stat_MeanNan(AV_biasfield(pGM));
    AV_biasfield = ScaleFactor .* AV_biasfield;

    % Put average images in RefIM
    % Reference image 4: average mean CBF
    RefIM{4} = repmat(xASL_im_rotate(AV_biasfield(:,:,53),90),[1 nSites]);





    %% ------------------------------------------------------------------------------------------------------------
    %% 4. Create new bias-fields, masked & expand
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




    %% ------------------------------------------------------------------------------------------------------------
    %% 5. Save biasfields
    % Store reference images
    RefIM{5} = xASL_im_rotate(BiasFieldMultipl(:,:,53,1),90);
    for iSite=2:nSites
        % Reference image 5: Biasfields
        RefIM{5} = [RefIM{5} xASL_im_rotate(BiasFieldMultipl(:,:,53,iSite),90)];
    end


    % Save biasfields
    for iSite=1:nSites
        FieldFileName = fullfile(x.D.PopDir,['Biasfield_Multipl_Site_' num2str(iSite) '_' x.S.SetsOptions{iSiteSet}{iSite} '.nii']);
        xASL_io_SaveNifti(pGMname,FieldFileName,BiasFieldMultipl(:,:,:,iSite));
    end



    %% ------------------------------------------------------------------------------------------------------------
    %% 6. Create backup of each CBF image before applying biasfields
    TryN = 1;
    BackupDir = fullfile(x.D.PopDir, 'BackupBeforeSiteRescale_1');
    while exist(BackupDir, 'dir')
          TryN = TryN+1;
          BackupDir = fullfile(x.D.PopDir, ['BackupBeforeSiteRescale_' num2str(TryN)]);
    end
    xASL_adm_CreateDir(BackupDir);
    fprintf('%s','Backing up CBF maps:   ');
    for iSub=1:x.dataset.nSubjects
        xASL_TrackProgress(iSub,x.dataset.nSubjects);
        for iSess=1:x.dataset.nSessions

            FilePath = fullfile(x.D.PopDir,[CBF_prefix '_' x.SUBJECTS{iSub} '_' x.SESSIONS{iSess} '.nii']);
            if xASL_exist(FilePath,'file')
                [~, File, Ext] = xASL_fileparts(FilePath);
                NewPath = fullfile(BackupDir,[File Ext]);
                xASL_Copy(FilePath, NewPath);
            end
        end
    end


    %% ------------------------------------------------------------------------------------------------------------
    %% 7. Print QC reference images
    xASL_vis_Imwrite([RefIM{1};RefIM{2};RefIM{3};RefIM{4}], fullfile(BackupDir,'Overview_biasfields.jpg') );

    
    %% ------------------------------------------------------------------------------------------------------------
    %% 8. Rescale CBF images (apply the biasfield correction)
    for iSite=1:nSites
        clear SiteScans nScans
        SiteScans = find(x.S.SetsID(:,iSiteSet)==AllSites(iSite));
        nScans = length(SiteScans);

        fprintf('\n%s',['Rescaling CBF images for intensity biasfield, site ' SiteNames{iSite} ':   ']);
        for iScan=1:nScans
            xASL_TrackProgress(iScan,nScans);
            clear iSubSess iSub iSess FileName tNII tIM
            iSubSess = SiteScans(iScan);
            iSub = ceil(iSubSess/x.dataset.nSessions);
            iSess = iSubSess- ((iSub-1)*x.dataset.nSessions);

            FileName = fullfile(x.D.PopDir,[CBF_prefix '_' x.SUBJECTS{iSub} '_' x.SESSIONS{iSess} '.nii']);
            if xASL_exist(FileName,'file')
                xASL_io_SaveNifti(FileName,FileName,xASL_io_Nifti2Im(FileName).*BiasFieldMultipl(:,:,:,iSite),[],0); % save rescaled image
            end
        end
    end

%%
%% Now normalize the SD for the sites (i.e. spatial CoV)




%% ------------------------------------------------------------------------------------------------------------
%% 9. Re-load & store site-average CBF images
IM = cell.empty(0);

iNext = ones(1,nSites);
for iSite=1:nSites
    SiteScans = find(x.S.SetsID(:,iSiteSet)==AllSites(iSite));
    nScans = length(SiteScans);
    fprintf('\n%s',['Loading CBF images for site ' SiteNames{iSite} ':   ']);
    for iScan=1:nScans
        xASL_TrackProgress(iScan,nScans);
        iSubSess = SiteScans(iScan);
        iSub = ceil(iSubSess/x.dataset.nSessions);
        iSess = iSubSess- ((iSub-1)*x.dataset.nSessions);
        FileName = fullfile(x.D.PopDir,[CBF_prefix '_' x.SUBJECTS{iSub} '_' x.SESSIONS{iSess} '.nii']);
        if xASL_exist(FileName,'file')
            IM{iSite}(:,iNext(iSite)) = xASL_im_IM2Column(xASL_io_Nifti2Im(FileName),x.S.masks.WBmask);
            NameList{iNext(iSite),iSite} = [x.SUBJECTS{iSub} '_' x.SESSIONS{iSess}];
            iNext(iSite) = iNext(iSite)+1;
        end
    end
end



%% ------------------------------------------------------------------------------------------------------------
%% 10. Perform site-wise SD correction & save new CBF images
% Divide by mean
for iSite=1:nSites
    MeanCBF{iSite} = squeeze(xASL_stat_MeanNan(IM{iSite},1)); % compute mean CBF
    IM{iSite} = IM{iSite}./MeanCBF{iSite}; % divide by mean CBF
    IM{iSite} = IM{iSite}-1; % demeaning
    StdCBF{iSite} = squeeze(xASL_stat_StdNan(IM{iSite},[],1));
    MeanStdCBF(iSite) = mean(StdCBF{iSite}(:));
end

MeanSiteStdRatio = MeanStdCBF./mean(MeanStdCBF); % Get average SD per scanner


% Correct images
for iSite=1:nSites
    IM{iSite} = IM{iSite}./MeanSiteStdRatio(iSite); % correct for average SD per scanner
    IM{iSite} = IM{iSite}+1; % undo demeaning
    IM{iSite} = IM{iSite}.*MeanCBF{iSite}; % multiply by mean CBF
end

% Save images
for iSite=1:nSites
    fprintf('\n%s',['Saving CBF images for site ' SiteNames{iSite} ':   ']);
    for iScan=1:size(IM{iSite},2)
        xASL_TrackProgress(iScan,nScans);
        FilePath = fullfile(x.D.PopDir,[CBF_prefix '_' NameList{iScan,iSite} '.nii']);
        xASL_io_SaveNifti(FilePath,FilePath,xASL_im_Column2IM(IM{iSite}(:,iScan),x.S.masks.WBmask),[],0);
    end
end


end