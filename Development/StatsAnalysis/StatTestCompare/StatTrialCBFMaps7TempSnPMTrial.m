%% Newer version, with slightly improved signal creation & ROI visualization.
%% GM masking implemented

%% CREATE two groups of maps, with 2nd group having BLOB
%   First fixed BLOB of 150% signal intensity
%   Using fixed variance among subjects (later take ROI variance)









%% Administration
x.StatsCmpDir    = fullfile( x.D.PopDir , 'StatsCompare');
xASL_adm_CreateDir(x.StatsCmpDir);



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
iInt                = 1 *50;
iVol                = 1 *50;
iSub                = 15;
iMask               = [14]; % [6 8 10 14];

% WhichIt             = [1 0 0];
% NoItSetting         = [75 33 15];   % settings to use when not iterating
% YesItMin            = [1 1 1];      % minimum iteration nr (start)
% YesItMax            = [100 100 15]; % maximum iteration nr (end)

x.StatsCmpDir     = fullfile( x.StatsCmpDir,'SnPM');
xASL_adm_CreateDir(x.StatsCmpDir);

tic
for iI=iInt  % intensities (effect height/size)
    Inten   =1+(iI/100);

    %% Create BLOB
    for iV=iVol  % volume size
        ActivationVolume                = iV/100; % percentage volume that is activated

        for iM=iMask

            %% Load images
            clear CBFImage
            for iIm=1:15 % groups
                CBFImage{1}(:,:,:,iIm)  = single(ASL.Data.data(PseudoRandom(iIm  ,1)*2,:,:,:)); % *2 = because of concatenated SubjectSessions (there are 2 sessions)
                CBFImage{2}(:,:,:,iIm)  = single(ASL.Data.data(PseudoRandom(iIm*2,1)*2,:,:,:));
            end

            [BLOB CBF_BLOB CBF1 GoldPos]    = Create3DROIblobOnData( x, CBFImage{2}, CBFImage{1}, iM, Inten, 0.3, ActivationVolume); % 12 = putamen. % cave, ASL input needs to be exact the same images!
            % Create3DROIblob uses the data to calculate the multiplication matrix needed to have a mean increase

            % Visualize
            % BLOB(isnan(BLOB))           = 0;
            % CBF_BLOB(isnan(CBF_BLOB))   = 0;
            % dip_image(xASL_im_rotate([background_mask(:,:,44)+BLOB(:,:,44);background_mask(:,:,44)+(CBFImage{1}(:,:,44,1)./50);background_mask(:,:,44)+(CBF_BLOB(:,:,44,1)./50)],90))
            % BLOB(BLOB==0)               = NaN;
            % CBF_BLOB(CBF_BLOB==0)       = NaN;
            % Plot2Dsurface(cropNaN(BLOB(:,:,44),0.1));
            % Plot2Dsurface(cropNaN(CBF_BLOB(:,:,44,1),0.1));
            % Plot2Dsurface(cropNaN(CBF_BLOB(:,:,44,2),0.1));
            % Plot2Dsurface(cropNaN(CBF_BLOB(:,:,44,3),0.1));


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


            % % Visualize
            % dip_image([xASL_im_rotate([CBFImage{1}(:,:,44,1);CBFImage{1}(:,:,44,2);CBFImage{1}(:,:,44,3)],90);xASL_im_rotate([CBFImage{2}(:,:,44,1);CBFImage{2}(:,:,44,2);CBFImage{2}(:,:,44,3)],90)])

            % diffIm                          = abs(CBFImage{1}-CBFImage{2});
            % dip_image(xASL_im_rotate([diffIm(:,:,44,1);diffIm(:,:,44,2);diffIm(:,:,44,3)],90))
            % dip_image(xASL_im_rotate(mean(diffIm(:,:,:,:),4),90))

            % xASL_im_CreateVisualFig_INPUT( mean(diffIm(:,:,:,:),4) )

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

                %% SnPM
                for CovMean=1:2 % covarying with mean no/yes
                    % First run
                    SnPMscript( x, iS       , CovMean-1, 10000, 1 ); % covarying with mean no/yes

                    % Where to store information
                    I(2)        = find(iI==iInt);
                    I(3)        = find(iV==iVol);
                    I(4)        = find(iM==iMask);
                    I(5)        = find(iS==iSub);

                    % 1) Uncorrected p=0.001 (unacceptable maybe)
                    % 2nd: get results (write filter image)
                    SnPMscript( x, iS       , CovMean-1, 10000, 0, 'uncorrected' );
                    StatResult                    = xASL_io_ReadNifti( fullfile( x.StatsCmpDir, 'SnPM_filtered.nii'));
                    StatResult                    = StatResult.dat(:,:,:); % Load SPM T-map
                    H                             = StatResult>0;

                    I(1)        = (CovMean*6)-5;
                    [TPR(I(1),I(2),I(3),I(4),I(5)) FPR(I(1),I(2),I(3),I(4),I(5))]   = CalcROC( GoldPos, H);

                    % 2) FWE p=0.05 based on peak-level (very strict)
                    SnPMscript( x, iS       , CovMean-1, 10000, 0, 'FWE' );
                    StatResult                    = xASL_io_ReadNifti( fullfile( x.StatsCmpDir, 'SnPM_filtered.nii'));
                    StatResult                    = StatResult.dat(:,:,:); % Load SPM T-map
                    H                             = StatResult>0;

                    I(1)        = (CovMean*6)-4;
                    [TPR(I(1),I(2),I(3),I(4),I(5)) FPR(I(1),I(2),I(3),I(4),I(5))]   = CalcROC( GoldPos, H);

                    % 3) FWE p=0.05 based on cluster-size (p=0.01 initially, FSL standard)
                    SnPMscript( x, iS       , CovMean-1, 10000, 0, 'cluster1' );
                    StatResult                    = xASL_io_ReadNifti( fullfile( x.StatsCmpDir, 'SnPM_filtered.nii'));
                    StatResult                    = StatResult.dat(:,:,:); % Load SPM T-map
                    H                             = StatResult>0;

                    I(1)        = (CovMean*6)-3;
                    [TPR(I(1),I(2),I(3),I(4),I(5)) FPR(I(1),I(2),I(3),I(4),I(5))]   = CalcROC( GoldPos, H);

                    % 4) FWE p=0.05 based on cluster-size (p=0.001 initially, SPM standard)
                    SnPMscript( x, iS       , CovMean-1, 10000, 0, 'cluster2' );
                    StatResult                    = xASL_io_ReadNifti( fullfile( x.StatsCmpDir, 'SnPM_filtered.nii'));
                    StatResult                    = StatResult.dat(:,:,:); % Load SPM T-map
                    H                             = StatResult>0;

                    I(1)        = (CovMean*6)-2;
                    [TPR(I(1),I(2),I(3),I(4),I(5)) FPR(I(1),I(2),I(3),I(4),I(5))]   = CalcROC( GoldPos, H);
                end

    %             AUC=TPR.*(1-FPR)

            end % subject
        end % mask
    end % volume
end % intensity
toc




% AUC(1)  = NaN;
figure(1);
% plot([1:length(AUC)]+100,AUC);
plot([1:length(AUC)],AUC(:,1),'b',[1:length(AUC)],AUC(:,2),'r');
% xlabel('CBF increase relatively to baseline (%)');
xlabel('Volume with respect to total ROI (%)');
ylabel('Area under ROC-curve (%)');
