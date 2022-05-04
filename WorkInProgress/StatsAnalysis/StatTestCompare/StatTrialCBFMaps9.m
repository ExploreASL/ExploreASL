%% Newer version, with slightly improved signal creation & ROI visualization.
%% GM masking implemented

%% CREATE two groups of maps, with 2nd group having BLOB
%   First fixed BLOB of 150% signal intensity
%   Using fixed variance among subjects (later take ROI variance)









%% Administration
x.StatsCmpDir    = fullfile( x.D.PopDir , 'StatsCompare');
xASL_adm_CreateDir(x.StatsCmpDir);

ROIs    = {'GM'	'WM'	'L-ICA'	'R-ICA'	'POS'	'Caudate'	'Cerebellum'	'Frontal'	'Insula'	'Occipital'	'Parietal' 'Putamen'	'Temporal' 'Thalamus'};

%% Pseudo-randomly CBF-maps 2nd session & divide into 2 groups
PseudoRandomFile    = fullfile( x.D.PopDir, 'PseudoRandomKey.mat');
if ~exist(PseudoRandomFile,'file')
    PseudoRandom(:,1)   = [1:30]';
    PseudoRandom(:,2)   = rand(1,30);
    PseudoRandom        = sortrows(PseudoRandom,2);
    save(PseudoRandomFile, 'PseudoRandom');
else load(PseudoRandomFile);
end


%% Load background for illustrating, & mask for masking

% GM	WM	L-ICA	R-ICA	POS	Caudate	Cerebellum	Frontal	Insula	Occipital	Parietal	Putamen	Temporal	Thalamus

background_mask     = fullfile( x.D.PopDir, 'DARTEL_mean_T1_template.nii');
background_mask     = xASL_io_ReadNifti( background_mask );
background_mask     = single( background_mask.dat(:,:,:) );
background_mask     = background_mask ./ max( background_mask(:) );

meanGMmask          = fullfile( x.D.PopDir, 'DARTEL_T1_template.nii');
meanGMmask          = nifti ( meanGMmask );
meanGMmask          = single( meanGMmask.dat(:,:,:) );
meanGMmask          = meanGMmask>(0.5*(max(meanGMmask(:))));
meanGMmask          = repmat( meanGMmask, [1 1 1 15]); % Cave 15 is subjectn!

%% Smoothing settings
% Fast matlab smoothing, visually slightly smaller than 8 mm FWHM in comparison with SPM smoothing, but nearly similar. KISS!!

FwHm                = 8/1.5;
FwHm2SD             = (2*(2*reallog(2))^0.5);
SD                  = FwHm/FwHm2SD;

%% To save iteration over 1) intensity, 2) volume, 3) masks, 4) nsubjects
% clear TPR FPR
nIt                 = [1 1 1 1];
iInt                = [20:40:100];
iVol                = [20:40:100];
iSub                = 15;
iMask               = [1 6]; % 10 14]; % occipital (large, cortical) & thalamus (small, subcortical)
SPMrun              = 0;
SnPMrun             = 1;

% ROI = GM	WM	L-ICA	R-ICA	POS	6 Caudate	Cerebellum	8 Frontal	Insula	10 Occipital	Parietal Putamen	Temporal 14	Thalamus

% WhichIt             = [1 0 0];
% NoItSetting         = [75 33 15];   % settings to use when not iterating
% YesItMin            = [1 1 1];      % minimum iteration nr (start)
% YesItMax            = [100 100 15]; % maximum iteration nr (end)


x.ResultDir   = fullfile( x.StatsCmpDir, 'Results');
xASL_adm_CreateDir( x.ResultDir );

tic

for iM=iMask
    for iI=iInt  % intensities (effect height/size)
        Inten   =1+(iI/100);
        for iV=iVol  % volume size
            ActivationVolume                = iV/100; % percentage volume that is activated

            %% Load images
            clear CBFImage BLOB CBF_BLOB CBF1 GoldPos
            for iIm=1:15 % groups
                CBFImage{1}(:,:,:,iIm)  = single(ASL.Data.data(PseudoRandom(iIm  ,1)*2,:,:,:)); % *2 = because of concatenated SubjectSessions (there are 2 sessions)
                CBFImage{2}(:,:,:,iIm)  = single(ASL.Data.data(PseudoRandom(iIm*2,1)*2,:,:,:));
            end

            % Create BLOB
            [BLOB CBF_BLOB CBF1 GoldPos]    = Create3DROIblobOnData( x, CBFImage{2}, CBFImage{1}, iM, Inten, 0.3, ActivationVolume); % 12 = putamen. % cave, ASL input needs to be exact the same images!
            % Create3DROIblob uses the data to calculate the multiplication matrix needed to have a mean increase

            %% Save descriptives PRE-BLOB (for all subjects)
            MaskedGoldPos                   = GoldPos & meanGMmask(:,:,:,1); % GoldPos is not masked yet, this happens after the smoothing below. But for comparison of "Outside_preBLOB" & "Outside_postBLOB", this MaskedGoldpos is required
            Diff                            = CBFImage{2}-CBFImage{1};
            MeanDiff                        = xASL_stat_MeanNan(Diff,4);
            StdDiff                         = xASL_stat_StdNan(Diff,0,4); % map of std over subjects (simplifying to paired comparison)

            temp                            = MeanDiff(MaskedGoldPos);
            desc_stats.Diff_GoldPos_preBLOB            = xASL_stat_MeanNan(temp(:));
            desc_stats.Stdspat_GoldPos_preBLOB         = xASL_stat_StdNan(temp(:)); % spat = spatial std on mean difference image
            temp                            = StdDiff(MaskedGoldPos);
            desc_stats.StdBS_GoldPos_preBLOB           = xASL_stat_MeanNan(temp(:));% BS = between subject

            OutsideGoldPosWithinGM          = ~MaskedGoldPos & meanGMmask(:,:,:,1);
            temp                            = MeanDiff(OutsideGoldPosWithinGM);
            desc_stats.Diff_Outside_preBLOB            = xASL_stat_MeanNan(temp(:));
            desc_stats.Stdspat_Outside_preBLOB         = xASL_stat_StdNan(temp(:)); % spat = spatial std on mean difference image
            temp                            = StdDiff(OutsideGoldPosWithinGM);
            desc_stats.StdBS_Outside_preBLOB           = xASL_stat_MeanNan(temp(:));% BS = between subject
            clear StdDiff MeanDiff Diff

            % Visualize
            % BLOB(isnan(BLOB))           = 0;
            % CBF_BLOB(isnan(CBF_BLOB))   = 0;
            % dip_image(xASL_im_rotate([background_mask(:,:,44)+BLOB(:,:,44);background_mask(:,:,44)+(CBFImage{1}(:,:,44,1)./50);background_mask(:,:,44)+(CBF_BLOB(:,:,44,1)./50)],90))
            % BLOB(BLOB==0)               = NaN;
            % CBF_BLOB(CBF_BLOB==0)       = NaN;
            % Plot2Dsurface(CropNaN(BLOB(:,:,44),0.1));
            % Plot2Dsurface(CropNaN(CBF_BLOB(:,:,44,1),0.1));
            % Plot2Dsurface(CropNaN(CBF_BLOB(:,:,44,2),0.1));
            % Plot2Dsurface(CropNaN(CBF_BLOB(:,:,44,3),0.1));


            %% Apply blob
            MaskTemp                        = repmat(GoldPos,[1 1 1 size(CBF_BLOB,4)]);
            CBFImage{2}(MaskTemp)           = CBF_BLOB((MaskTemp));
            clear MaskTemp


            %% Smooth
            for iIm     = 1:size(CBFImage{1},4)
                %CBFImage{1}(:,:,:,iIm)      = dip_array(gaussf ( CBFImage{1}(:,:,:,iIm) , [SD SD SD], 'fir'));
                %CBFImage{2}(:,:,:,iIm)      = dip_array(gaussf ( CBFImage{2}(:,:,:,iIm) , [SD SD SD], 'fir'));
				CBFImage{1}(:,:,:,iIm)      = xASL_im_ndnanfilter(CBFImage{1}(:,:,:,iIm),'gauss',[SD SD SD]*2.335,0);
				CBFImage{2}(:,:,:,iIm)      = xASL_im_ndnanfilter(CBFImage{2}(:,:,:,iIm),'gauss',[SD SD SD]*2.335,0);
            end

            %   Visualize
            %     close all
            %     dip_image(BLOB(:,:,51))
            %     dip_image(CBF_BLOB(:,:,51))
            %     dip_image( CBFImage{2}(:,:,51,1) )
            %     dip_image( mean( CBFImage{2}(:,:,51,:) - CBFImage{1}(:,:,51,:) ,4 ) )


            %% Mask GM
            CBFImage{1}(~meanGMmask)        = NaN;
            CBFImage{2}(~meanGMmask)        = NaN;
            GoldPos(~meanGMmask(:,:,:,1))   = 0;

            %% Save descriptives POST-BLOB (for all subjects)
            % Get representative slice for ROI
            xASL_vis_CreateVisualFig_INPUT( [GoldPos(:,:,:,1)], fullfile(x.ResultDir,['GoldPos_M' num2str(iM) '_I' num2str(iI) '_V' num2str(iV)  '.jpg']) );
            xASL_vis_CreateVisualFig_INPUT( BLOB(:,:,:,1), fullfile(x.ResultDir,['BLOB_M' num2str(iM) '_I' num2str(iI) '_V' num2str(iV)  '.jpg']) );
            Diff        = CBFImage{2}-CBFImage{1};
            MeanDiff    = xASL_stat_MeanNan(Diff,4);
            StdDiff                         = xASL_stat_StdNan(Diff,0,4); % map of std over subjects (simplifying to paired comparison)

            xASL_vis_CreateVisualFig_INPUT( abs(MeanDiff), fullfile( x.ResultDir,['MeanDiff_M' num2str(iM) '_I' num2str(iI) '_V' num2str(iV)  '.jpg']) );
            xASL_vis_CreateVisualFig_INPUT( StdDiff , fullfile( x.ResultDir,['StdDiff_M' num2str(iM) '_I' num2str(iI) '_V' num2str(iV)  '.jpg']) );

            % Save descriptives
            temp                            = MeanDiff(MaskedGoldPos);
            desc_stats.Diff_GoldPos_postBLOB           = xASL_stat_MeanNan(temp(:));
            desc_stats.Stdspat_GoldPos_postBLOB        = xASL_stat_StdNan(temp(:)); % spat = spatial std on mean difference image
            temp                            = StdDiff(MaskedGoldPos);
            desc_stats.StdBS_GoldPos_postBLOB          = xASL_stat_MeanNan(temp(:));% BS = between subject

            temp                            = MeanDiff(OutsideGoldPosWithinGM);
            desc_stats.Diff_Outside_postBLOB           = xASL_stat_MeanNan(temp(:));
            desc_stats.Stdspat_Outside_postBLOB        = xASL_stat_StdNan(temp(:)); % spat = spatial std on mean difference image
            temp                            = StdDiff(OutsideGoldPosWithinGM);
            desc_stats.StdBS_Outside_postBLOB          = xASL_stat_MeanNan(temp(:));% BS = between subject

            SaveFile                        = fullfile(x.ResultDir, ['Desc_stats_M' num2str(iM) '_I' num2str(iI) '_V' num2str(iV)  '_S15.mat']);
            save(SaveFile,'desc_stats');
            clear SaveFile desc_stats temp
            clear StdDiff MeanDiff Diff

            % % Visualize
            % dip_image([xASL_im_rotate([CBFImage{1}(:,:,44,1);CBFImage{1}(:,:,44,2);CBFImage{1}(:,:,44,3)],90);xASL_im_rotate([CBFImage{2}(:,:,44,1);CBFImage{2}(:,:,44,2);CBFImage{2}(:,:,44,3)],90)])

            % diffIm                          = abs(CBFImage{1}-CBFImage{2});
            % dip_image(xASL_im_rotate([diffIm(:,:,44,1);diffIm(:,:,44,2);diffIm(:,:,44,3)],90))
            % dip_image(xASL_im_rotate(mean(diffIm(:,:,:,:),4),90))

            % xASL_vis_CreateVisualFig_INPUT( mean(diffIm(:,:,:,:),4) )

            %% Save
            % First delete whole directory
            xASL_adm_DeleteFileList(x.StatsCmpDir, '^.*\..*$', 0);

            ExampleNii      = fullfile( x.D.PopDir, ['DARTEL_c1T1_' x.SUBJECTS{1} '.nii']);

            for iIm     = 1:size(CBFImage{1},4)
                NewName                 = fullfile( x.StatsCmpDir, ['CBF1_' num2str(iIm) '.nii']);
                xASL_io_SaveNifti( ExampleNii, NewName, CBFImage{1}(:,:,:,iIm) );
                clear NewName
                NewName     = fullfile( x.StatsCmpDir, ['CBF2_' num2str(iIm) '.nii']);
                xASL_io_SaveNifti( ExampleNii, NewName, CBFImage{2}(:,:,:,iIm) );
                clear NewName
            end

            % xASL_stat_SumNan(xASL_stat_SumNan(xASL_stat_SumNan( CBFImage{2}(:,:,:,1)>0 ))) voxel count is correct

            %% Iterating number of subjects (sample size)
            nSubjects   =15;

            for iS=iSub % subjects to include
    %             iSub=nSubjects-iSub+1;
                % Remove subject maps

                % Check if exists & if requested
                SaveFile                        = fullfile(x.ResultDir, ['SPM_M' num2str(iM) '_I' num2str(iI) '_V' num2str(iV)  '_S' num2str(iS) '.mat']);

                if ~exist(SaveFile,'file') && SPMrun

                    %% Standard SPM
                    for CovMean=1:2 % covarying with mean no/yes

                        SPMDesignT2Loop( x, iS, CovMean-1 ); % covarying with mean no/yes

                        StatResult                      = xASL_io_ReadNifti( fullfile( x.StatsCmpDir, 'spmT_0001.nii'));
                        StatResult                      = StatResult.dat(:,:,:); % Load SPM T-map

    %                     % Where to store information
    %                     I(2)        = find(iI==iInt);
    %                     I(3)        = find(iV==iVol);
    %                     I(4)        = find(iM==iMask);
    %                     I(5)        = find(iS==iSub);

                        % 1) Uncorrected p=0.001 (unacceptable maybe)
                        clear StatResult H
                        H           = StatResult>3.4;

                        I(1)        = (CovMean*4)-3;
                        [TPR(I(1)) FPR(I(1))]           = CalcROC( GoldPos, H);

                        % 2) FWE p=0.05 based on peak-level (very strict)
                        clear StatResult H
                        H           = StatResult>5.7;

                        I(1)        = (CovMean*4)-2;
                        [TPR(I(1)) FPR(I(1))]           = CalcROC( GoldPos, H);

                        % 3) FWE p=0.05 based on cluster-size (p=0.01 initially, FSL standard)
                        clear StatResult H
                        SPMCorrClusTh( x, 0.01, 0.05 );
                        StatResult                    = xASL_io_ReadNifti( fullfile( x.StatsCmpDir, 'FilterSPM.nii'));
                        StatResult                    = StatResult.dat(:,:,:); % Load SPM T-map
                        H                             = StatResult>0;

                        I(1)        = (CovMean*4)-1;
                        [TPR(I(1)) FPR(I(1))]           = CalcROC( GoldPos, H);

                        % 4) FWE p=0.05 based on cluster-size (p=0.001 initially, SPM standard)
                        clear StatResult H
                        SPMCorrClusTh( x, 0.001, 0.05 );
                        StatResult                    = xASL_io_ReadNifti( fullfile( x.StatsCmpDir, 'FilterSPM.nii'));
                        StatResult                    = StatResult.dat(:,:,:); % Load SPM T-map
                        H                             = StatResult>0;

                        I(1)        = (CovMean*4)-0;
                        [TPR(I(1)) FPR(I(1))]           = CalcROC( GoldPos, H);
                        % 5) FWE p=0.05 based on cluster-size defined by TFCE (threshold-free cluster estimation)
        %                 Index(1)      = (CovMean*5)-0;
            %             TFCscript( x, 500 );
            %             StatResult                  = xASL_io_ReadNifti( fullfile( x.StatsCmpDir, 'FilterSPM.nii'));
            %             StatResult                  = StatResult.dat(:,:,:); % Load SPM T-map
            %             [TPR(ID) FPR(ID)]             = CalcROC( GoldPos, StatResult>0);


                        % Delete temp files
                        xASL_adm_DeleteFileList(x.StatsCmpDir,'^beta_.*$');
                        delete('beta_0001.nii');
                        delete('beta_0002.nii');
                        delete('con_0001.nii');
                        delete('mask.nii');
                        delete('ResMS.nii');
                        delete('RPV.nii');
                        delete('SPM.mat');
                        delete('spmT_0001.nii');
                        xASL_adm_DeleteFileList(x.StatsCmpDir,'^FilterSPM\.(nii|nii\.gz)$');
                    end

                    % Save data
                    save( SaveFile, 'TPR','FPR' );
                    clear TPR FPR
                end


                % Check if exists & if requested
                SaveFile                        = fullfile(x.ResultDir, ['SnPM_M' num2str(iM) '_I' num2str(iI) '_V' num2str(iV)  '_S' num2str(iS) '.mat']);
                if ~exist(SaveFile,'file') && SnPMrun

                    %% SnPM
                    for CovMean=1:2 % covarying with mean no/yes

                        % First run
                        SnPMscript( x, iS       , CovMean-1, 10000, 1 ); % covarying with mean no/yes

                        % 1) Uncorrected p=0.001 (unacceptable maybe)
                        % 2nd: get results (write filter image)
                        clear StatResult H
                        SnPMscript( x, iS       , CovMean-1, 1000, 0, 'uncorrected' );
                        StatResult                    = xASL_io_ReadNifti( fullfile( x.StatsCmpDir, 'SnPM_filtered.nii'));
                        StatResult                    = StatResult.dat(:,:,:); % Load SPM T-map
                        H                             = StatResult>0;

                        I(1)        = (CovMean*4)-3;
                        [TPR(I(1)) FPR(I(1))]           = CalcROC( GoldPos, H);

                        % 2) FWE p=0.05 based on peak-level (very strict)
                        clear StatResult H
                        SnPMscript( x, iS       , CovMean-1, 10000, 0, 'FWE' );
                        StatResult                    = xASL_io_ReadNifti( fullfile( x.StatsCmpDir, 'SnPM_filtered.nii'));
                        StatResult                    = StatResult.dat(:,:,:); % Load SPM T-map
                        H                             = StatResult>0;

                        I(1)        = (CovMean*4)-2;
                        [TPR(I(1)) FPR(I(1))]           = CalcROC( GoldPos, H);

                        % 3) FWE p=0.05 based on cluster-size (p=0.01 initially, FSL standard)
                        clear StatResult H
                        SnPMscript( x, iS       , CovMean-1, 10000, 0, 'cluster1' );
                        StatResult                    = xASL_io_ReadNifti( fullfile( x.StatsCmpDir, 'SnPM_filtered.nii'));
                        StatResult                    = StatResult.dat(:,:,:); % Load SPM T-map
                        H                             = StatResult>0;

                        I(1)        = (CovMean*4)-1;
                        [TPR(I(1)) FPR(I(1))]           = CalcROC( GoldPos, H);

                        % 4) FWE p=0.05 based on cluster-size (p=0.001 initially, SPM standard)
                        clear StatResult H
                        SnPMscript( x, iS       , CovMean-1, 10000, 0, 'cluster2' );
                        StatResult                    = xASL_io_ReadNifti( fullfile( x.StatsCmpDir, 'SnPM_filtered.nii'));
                        StatResult                    = StatResult.dat(:,:,:); % Load SPM T-map
                        H                             = StatResult>0;

                        I(1)        = (CovMean*4)-0;
                        [TPR(I(1)) FPR(I(1))]           = CalcROC( GoldPos, H);

                        % Delete temp files
                        xASL_adm_DeleteFileList(x.StatsCmpDir,'^beta_.*$');
                        xASL_adm_DeleteFileList(x.StatsCmpDir,'^lP.*$');
                        xASL_adm_DeleteFileList(x.StatsCmpDir,'^ResMS.*$');
                        xASL_adm_DeleteFileList(x.StatsCmpDir,'^SnPM.*$');
                        xASL_adm_DeleteFileList(x.StatsCmpDir,'^snpm.*$');
                        xASL_adm_DeleteFileList(x.StatsCmpDir,'^spm.*$');
                        xASL_adm_DeleteFileList(x.StatsCmpDir,'^XYZ.*$');
                        xASL_adm_DeleteFileList(x.StatsCmpDir,'^SnPM_filtered.*$');
                        xASL_adm_DeleteFileList(x.StatsCmpDir,'^STCS.*$');
                    end

                    % Save data
                    save( SaveFile, 'TPR','FPR' );
                    clear TPR FPR
                end

            end
        end
    end
end
toc

%% CAVE change permutations SnPM!!!


%%
%%
%%
%%
%%

%% Load and write descriptive statistics





%% Load & show results
clear iTech iInt iMask iVol iSub iT iI iV iS iM iTread iIread iVread iSread iMread
for iInt=[20:20:100]
    for iMask=[10 14] % Cave: these are iterations of the plots that we want to show! When "for" is used, the total vector is not used but rather "looped over"

        PlotFactor          = 'iVol';
        iVol                = [20:20:100]; % These are iterations within a plot!
        % iInt                = 20;
        iSub                = 15;
    %     iMask               = 10; % occipital (large, cortical) & thalamus (small, subcortical)
        iTech               = {'SPM'} ;% SnPM}';

        % Create file for descriptive stats
        switch PlotFactor
        case 'iInt'
            FileDescStats   = fullfile( x.ResultDir, ['Desc_stats_' ROIs{iMask} '_iterate_intensity_V' num2str(iVol) '_S' num2str(iSub) '.csv']);
        case 'iVol'
            FileDescStats   = fullfile( x.ResultDir, ['Desc_stats_' ROIs{iMask} '_iterate_volume_I'    num2str(iInt) '_S' num2str(iSub) '.csv']);
        case 'iSub' % not yet supported
            FileDescStats   = fullfile( x.ResultDir, ['Desc_stats_' ROIs{iMask} '_iterate_subjects_I'  num2str(iInt) '_V' num2str(iVol) '.csv']);
        end

        x.S.summary_fid        = fopen(FileDescStats,'wt');
        % Print header

        switch PlotFactor
        case 'iInt'
            fprintf(x.S.summary_fid,'%s,','Intensity (%)');
        case 'iVol'
            fprintf(x.S.summary_fid,'%s,','Volume (%)');
        case 'iSub' % not yet supported
            fprintf(x.S.summary_fid,'%s,','Sample size (n)');
        end

        fprintf(x.S.summary_fid,'%s','Diff_GoldPos_preBLOB   ,Stdspat_GoldPos_preBLOB  ,StdBS_GoldPos_preBLOB,');
        fprintf(x.S.summary_fid,'%s','Diff_Outside_preBLOB   ,Stdspat_Outside_preBLOB  ,StdBS_Outside_preBLOB,');
        fprintf(x.S.summary_fid,'%s','Diff_GoldPos_postBLOB  ,Stdspat_GoldPos_postBLOB ,StdBS_GoldPos_postBLOB,');
        fprintf(x.S.summary_fid,'%s','Diff_Outside_postBLOB  ,Stdspat_Outside_postBLOB ,StdBS_Outside_postBLOB,');
        fprintf(x.S.summary_fid,'\n');


        % Load results
        clear TPRtotal FPRtotal
        for iT=1:length(iTech)
            iTread=iTech{iT};
            for iM=1:length(iMask)
                iMread=iMask(iM);
                for iI=1:length(iInt)
                    iIread=iInt(iI);
                    for iV=1:length(iVol)
                        iVread=iVol(iV);
                        for iS=1:length(iSub)
                            iSread=iSub(iS);

                            % Create TPR & FPR overview
                            clear TPR FPR FileLoad
                            FileLoad        = fullfile( x.ResultDir, [iTread '_M' num2str(iMread) '_I' num2str(iIread) '_V' num2str(iVread) '_S' num2str(iSread) '.mat']);
                            load( FileLoad )
                            TPRtotal(:,iT,iM,iI,iV,iS)  = TPR';
                            FPRtotal(:,iT,iM,iI,iV,iS)  = FPR';

                            % Load & print descriptive stats
                            clear FileLoad desc_stats
                            FileLoad        = fullfile( x.ResultDir, ['Desc_stats_M' num2str(iMread) '_I' num2str(iIread) '_V' num2str(iVread) '_S' num2str(iSread) '.mat']);
                            load( FileLoad );

                            switch PlotFactor % first cell of row states iteration number
                            case 'iInt'
                                fprintf(x.S.summary_fid,'%s,',num2str(iIread));
                            case 'iVol'
                                fprintf(x.S.summary_fid,'%s,',num2str(iVread));
                            case 'iSub' % not yet supported
                                fprintf(x.S.summary_fid,'%s,',num2str(iSread));
                            end

                            % Print descriptive statistic values
                            fprintf(x.S.summary_fid,'%s,%s,%s,',num2str(desc_stats.Diff_GoldPos_preBLOB) ,num2str(desc_stats.Stdspat_GoldPos_preBLOB) ,num2str(desc_stats.StdBS_GoldPos_preBLOB ) );
                            fprintf(x.S.summary_fid,'%s,%s,%s,',num2str(desc_stats.Diff_Outside_preBLOB) ,num2str(desc_stats.Stdspat_Outside_preBLOB) ,num2str(desc_stats.StdBS_Outside_preBLOB ) );
                            fprintf(x.S.summary_fid,'%s,%s,%s,',num2str(desc_stats.Diff_GoldPos_postBLOB),num2str(desc_stats.Stdspat_GoldPos_postBLOB),num2str(desc_stats.StdBS_GoldPos_postBLOB) );
                            fprintf(x.S.summary_fid,'%s,%s,%s,',num2str(desc_stats.Diff_Outside_postBLOB),num2str(desc_stats.Stdspat_Outside_postBLOB),num2str(desc_stats.StdBS_Outside_postBLOB) );
                            fprintf(x.S.summary_fid,'\n');

                        end
                    end
                end
            end
        end

        % Close descriptive statistic file
        fclose(x.S.summary_fid);

        % Assemble the TPR & FPR
        TPRtotal        = squeeze(TPRtotal);
        FPRtotal        = squeeze(FPRtotal);

        close all
        plotColorOptions        = {'r'      'b'     'y'         'k'     'c'     'm'         'g'     'r' };
        plotColorNames          = {'red'    'blue'  'yellow'    'black' 'cyan'  'magenta'   'green'  'red'};

        switch PlotFactor
        case 'iInt'
            printX      = 'Effect size (%)';
            plotX       = iInt;
            PrintTitle  = ['Iteration_Intensity_' char(iTech) '_' ROIs{iMask} '_V' num2str(iVol) '_S' num2str(iSub)];
        case 'iVol'
            printX      = 'Volume with respect to total ROI (%)';
            plotX       = iVol;
            PrintTitle  = ['Iteration_Volume_' char(iTech) '_' ROIs{iMask} '_I' num2str(iInt) '_S' num2str(iSub)];
        case 'iSub'
            printX      = 'Sample size (n)';
            plotX       = iSub;
            PrintTitle  = ['Iteration_Subject_' char(iTech) '_' ROIs{iMask} '_I' num2str(iInt) '_V' num2str(iVol)];
        end

        fig = figure('Visible','off');
        subplot(2,1,1)
        hold on
        for iPlot=1:4 % no co-varying mean
            CO  = [plotColorOptions{iPlot} '-'];
            plot(plotX,(TPRtotal(iPlot,:) ),CO) % Volume
        end
        for iPlot=5:size(TPRtotal,1) % co-varying mean
            CO  = [plotColorOptions{iPlot-4} ':'];
            plot(plotX,(TPRtotal(iPlot,:) ),CO) % Volume
        end

        ylabel('True positive rate (%)');
        xlabel(printX);
        V=axis;
        axis([V(1:2) 0 1]);
        h   = title( PrintTitle  );
        set(h,'interpreter','none');
        subplot(2,1,2)
        hold on

        for iPlot=1:4 % no co-varying mean
            CO  = [plotColorOptions{iPlot} '-'];
            plot(plotX,squeeze(FPRtotal(iPlot,:) ),CO) % Volume
        end
        for iPlot=5:size(FPRtotal,1) % co-varying mean
            CO  = [plotColorOptions{iPlot-4} ':'];
            plot(plotX,squeeze(FPRtotal(iPlot,:) ),CO) % Volume
        end

        ylabel('False positive rate (%)');
        xlabel(printX);
        V=axis;
        axis([V(1:2) 0 0.1]);
        h   = title( PrintTitle  );
        set(h,'interpreter','none');
        SaveFile    = fullfile( x.ResultDir, PrintTitle );
        saveas( fig ,SaveFile,'jpg');
        clear SaveFile fig PrintTitle
        close


    end
end
