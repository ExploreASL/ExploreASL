function xASL_wrp_RealignASL(x,bSubtraction)
%xASL_wrp_RealignASL Submodule of ExploreASL ASL Module, that realigns
%volumes
%
% FORMAT: xASL_wrp_RealignASL(x[,bSubtraction])
%
% INPUT:
%   x               - structure containing fields with all information required to run this submodule (REQUIRED)
%   bSubtraction    - boolean for subtraction-based imaging (true, ASL) or not (false, e.g. fMRI, DTI etc). (OPTIONAL, DEFAULT = true) 
%
% OUTPUT: n/a (registration changes the NIfTI orientation header only,
%              with the exception of the affine transformation, which is
%              saved separately as x.P.Path_mean_PWI_Clipped_sn_mat
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This submodule estimates motion by spm_realign, which uses a
% rigid-body registration (3 translations, 3 rotations). It runs ENABLE to
% reject outliers and provides a visualization. ENABLE, QC and visualizations
% are based on the Net Displacement Vector (NDV) (in mm):
% according to Pythagorean/Euclydian RMS
% https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind1211&L=fsl&P=R34458&1=fsl&9=A&J=on&d=No+Match%3BMatch%3BMatches&z=4
% view this link for image of rotation roll, pitch and yaw https://www.google.nl/search?q=rotation+pitch+yaw+roll&espv=2&tbm=isch&imgil=LW3Nn1K-L6Oc7M%253A%253B-aSyykkRityJoM%253Bhttp%25253A%25252F%25252Fwww.grc.nasa.gov%25252FWWW%25252Fk-12%25252Fairplane%25252Frotations.html&source=iu&usg=__MlLQ5VuyRbm6kZP0vBJlPxmfbkw%3D&sa=X&ei=TWfjU4WcK4bqyQPqu4Fo&ved=0CD8Q9QEwBQ&biw=1680&bih=946#facrc=_&imgdii=_&imgrc=LW3Nn1K-L6Oc7M%253A%3B-aSyykkRityJoM%3Bhttp%253A%252F%252Fwww.grc.nasa.gov%252FWWW%252Fk-12%252Fairplane%252FImages%252Frotations.gif%3Bhttp%253A%252F%252Fwww.grc.nasa.gov%252FWWW%252Fk-12%252Fairplane%252Frotations.html%3B709%3B533
% 
% This submodule performs the following steps:
%
% 1. Estimate motion
% 2. Calculate and plot position and motion parameters
% 3. Threshold-free spike definition (based on ENABLE, but with t-stats rather than the threshold p<0.05)
% 4. Remove spike frames from nifti
%
% EXAMPLE: xASL_wrp_RealignASL(x);
% __________________________________
% Copyright (c) 2015-2021 ExploreASL


%% ----------------------------------------------------------------------------------------
%% Administration

if ~exist('bSubtraction','var')
    bSubtraction = 1; % this tells script that timeseries contain subtractive/pair-wise data
end

if  bSubtraction
    InputPath = x.P.Path_ASL4D;
else
    InputPath = x.P.Path_func_bold; % or DTI?
end

[Fpath, Ffile, Fext] = fileparts(InputPath);
rpfile = fullfile( Fpath, ['rp_' Ffile '.txt']);
rInputPath = fullfile( Fpath, ['r' Ffile Fext]);

if ~isfield(x.modules.asl,'SpikeRemovalThreshold')
    fprintf('%s\n','x.modules.asl.SpikeRemovalThreshold was not defined yet, default setting = 0.01 used');
    x.modules.asl.SpikeRemovalThreshold   = 0.01; % default threshold, decreased this from 0.05 to 0.01,
    % since we want to remove Spikes, perhaps except very small spikes
end

%% Set defaults
exclusion = NaN;
PercExcl = NaN;
MinimumtValue = NaN;

if x.modules.asl.bMultiPLD || x.modules.asl.bMultiTE || x.modules.asl.bDeltaM
    % ENABLE is disabled if multiPLD/TE
    bENABLE = 0;
    bZigZag = 0;
else
	if bSubtraction
		bENABLE = 1;
		bZigZag = 1;
	else
		bENABLE = 0;
		bZigZag = 0;
	end
end

%% ----------------------------------------------------------------------------------------
%% 1 Estimate motion
fprintf('SPM motion estimation');

% Define realignment settings
tempnii = xASL_io_ReadNifti(InputPath);
nFrames = double(tempnii.hdr.dim(5));
if length(x.EchoTime)>1
    nFrames=nFrames/numel(unique(x.EchoTime));
end
minVoxelSize = double(min(tempnii.hdr.pixdim(2:4)));

% Issue warning if empty image
if max(max(max(max(tempnii.dat(:)))))==0 || numel(unique(tempnii.dat(:)))==1
    warning('Invalid input image, skipping');
    return;
end

switch x.settings.Quality
	case 1 % normal quality
    flags.quality = 1;
    flags.sep = minVoxelSize;
    case 0 % low quality for fast try-out
    flags.quality = 0.01;
    flags.sep = minVoxelSize*2;
end

flags.rtm = 1; % realign to mean
flags.interp = 1;
flags.graphics = 0;

% If previous realign parameters exist, delete them
xASL_delete(rpfile);

% Run SPM

% Calculating the motion for every first TE of each PLD and
% considering it the same for the other TEs from that PLD.
%
% if nFrames>2 && bSubtraction && length(x.EchoTime)>1  %Multi TE 
%     uniqueTE=uniquetol(x.EchoTime); %gives the number of unique TEs
%     NumTEs=numel(uniqueTE);
%     minTE=min(uniqueTE); 
%     positionMinTE=find(x.EchoTime == minTE); %positions that have the min TE
%     ImInfo=spm_vol(InputPath);
%     ImInfoFirstTEs=ImInfo(positionMinTE);
%     
%     if numel(unique(x.Q.Initial_PLD))==1 %multiTE + single PLD
%         spm_realign(ImInfoFirstTEs,flags,true);
%         MotionFirstTEs=load(rpfile);
%         MotionAllTEs=repelem(MotionFirstTEs(:,:),NumTEs,1); %repeats each row NumTEs times
%         save('rp_ASL4D.txt','MotionAllTEs','-ascii') %saves the rpfile again into .txt
%         
%     elseif numel(unique(x.Q.Initial_PLD))>1 % multiTE + multiPLD=Hadamard
%         spm_realign(ImInfoFirstTEs,flags,false);
%         MotionFirstTEs=load(rpfile);
%         MotionAllTEs=repelem(MotionFirstTEs(:,:),NumTEs,1); 
%         save('rp_ASL4D.txt','MotionAllTEs','-ascii')
%     end

% Run motion correction for corresponding case
if nFrames>2 && bSubtraction && (x.modules.asl.bMultiPLD  ||  x.modules.asl.bMultiTE)
    % Multi-PLD, Multi-TE or Hadamard
    spm_realign(spm_vol(InputPath),flags,false);
    
elseif nFrames>2 && bSubtraction
    spm_realign(spm_vol(InputPath),flags,true);
    
elseif nFrames>1
    spm_realign(spm_vol(InputPath),flags,false);
end

%% ----------------------------------------------------------------------------------------
%% 2 Calculate and plot position and motion parameters
fprintf('%s\n','Calculate & plot position & motion parameters');

% Summarize real-world realign parameters into net displacement vector (NDV)
rp = load(rpfile, '-ascii'); % load the 3 translation and 3 rotation values
MeanRadius = 50; % typical distance center head to cerebral cortex (Power et al., NeuroImage 2012)
% PM: assess this from logical ASL EPI mask? This does influence the weighting of rotations compared to translations

% % FD = frame displacement
% if length(x.EchoTime)>2
%     FD{1} = rp(1:NumTEs:end,:); %multiTE -> gives back the normal rp for the plots
%     FD{2} = diff(rp(1:NumTEs:end,:));
% else
FD{1}=rp; % position (absolute displacement)
FD{2} = diff(rp); % motion (relative displacement)
% end


if max(rp(:))==0
    warning('Something wrong with motion parameters, skipping');
    return;
end

close all;
if usejava('jvm') % only if JVM loaded
    fig = figure('Visible','off');
else
    fprintf('Warning, skipping motion plot, JVM missing\n');
end

for ii=1:2
    tx{ii} = FD{ii}(:,1); ty{ii} = FD{ii}(:,2); tz{ii}  = FD{ii}(:,3); % translations
    rx{ii} = FD{ii}(:,4); ry{ii} = FD{ii}(:,5); rz{ii}  = FD{ii}(:,6); % rotations (pitch, roll, yaw)

    PartTranslation{ii} = tx{ii}.^2 + ty{ii}.^2 + tz{ii}.^2;
    PartRotation{ii} = 0.2*MeanRadius^2* ((cos(rx{ii})-1).^2 + (sin(rx{ii})).^2 + (cos(ry{ii})-1).^2 + (sin(ry{ii})).^2 + (cos(rz{ii})-1).^2 + (sin(rz{ii})).^2);
    try
        NDV{ii} = sqrt(PartTranslation{ii} + PartRotation{ii});
    catch
        
    end

    if ii==2
        NDV{2} = [0; NDV{2}]; % add leading zero difference
    end

    % Descriptives
    median_NDV{ii} = median(NDV{ii});
    mean_NDV{ii} = mean(NDV{ii});
    max_NDV{ii} = max(NDV{ii});
    SD_NDV{ii} = std(NDV{ii});
    MAD_NDV{ii} = xASL_stat_MadNan(NDV{ii},0); % median absolute deviation from median

    if usejava('jvm') % only if JVM loaded
        subplot(3,1,ii); % plot position (subplot 1) & motion (subplot 2)
        plot(NDV{ii},'Color',[0.4,0.4,0.4]); % lines between frames
        hold on
        plot(NDV{ii},'o','MarkerSize',5); % cirkels for frames
        hold on

        plot(repmat(mean_NDV{ii},151,1),'Color',[0,0,1]); % mean NDV in blue
        hold on

        if ii==1
            axis([0 nFrames 0 max(NDV{ii})]);
            title(['Position plot of ' x.P.SubjectID '-' x.P.SessionID ' relative to first frame']);
            ylabel('NDV (mm)');

        elseif ii==2

            axis([0 nFrames 0 minVoxelSize]);
            title(['Motion plot of ' x.P.SubjectID '-' x.P.SessionID]);
            ylabel('NDV/frame (mm//frame)');
        end

        xlabel('frame#');
        axis([1 length(NDV{ii}) 0 ceil(max(NDV{ii}))]); % fix X-axes to be same for subplots
    end
end

%% ----------------------------------------------------------------------------------------
%% 3) Threshold-free spike definition (based on ENABLE, but with t-stats rather than the threshold p<0.05)


if bSubtraction && nFrames<=10
    fprintf('Too few control-label pairs for ENABLE, skipping\n');
    bENABLE = 0;
end

if bENABLE && bSubtraction && nFrames>10 % == more than 5 pairs
    
    % Sort motion of control-label pairs
    MotionTime = NDV{2}; % motion
    MotionTime = MotionTime(1:2:end-1)+MotionTime(2:2:end); % additive motion for each control-label pair
    MotionTime(:,2) = [1:1:length(MotionTime)];
    MotionTimeSort = sortrows(MotionTime,1);
    
    % Resample ASL image (apply motion estimation)
    xASL_delete(rInputPath);
    matlabbatch{1}.spm.spatial.realign.write.data = {InputPath};
    matlabbatch{1}.spm.spatial.realign.write.roptions.which = [2 0];
    matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 1;
    matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'r';
    spm_jobman('run',matlabbatch);
    
    
    % Create a mask from the mean PWI
    xASL_io_PairwiseSubtraction(rInputPath, x.P.Path_mean_PWI_Clipped, 0, 0); % create PWI & mean_control
    MaskIm = xASL_im_ClipExtremes(x.P.Path_mean_PWI_Clipped,0.95,0.7);
    MaskIm = MaskIm>min(MaskIm(:));
    xASL_delete(x.P.Path_mean_PWI_Clipped);
    
    % Load ASL-image
    IM = xASL_io_Nifti2Im(rInputPath);
    
    if  bSubtraction % if subtractive/pairwise data
        [~, ~, OrderContLabl] = xASL_quant_GetControlLabelOrder(IM);
        
        if OrderContLabl~=1
            IM = IM(:,:,:,1:2:end-1)-IM(:,:,:,2:2:end); % control-label order doesn't matter for one-sample ttest p-values
        else
            IM = IM(:,:,:,2:2:end)-IM(:,:,:,1:2:end-1); % control-label order doesn't matter for one-sample ttest p-values
        end
    end
    
    
    fprintf('Running ENABLE:   ');
    tValue(1,1) = 0;
    SortIM = IM(:,:,:,MotionTimeSort(:,2)); % Sort ASL-pairs by motion
    for iVolume = 1:size(SortIM,4)
        TempIm = SortIM(:,:,:,iVolume);
        SortMask(:,iVolume) = TempIm(MaskIm); % create data columns
    end
    
    for iVolume=2:size(SortIM,4)
        xASL_TrackProgress(iVolume,size(SortMask,2));
        TempTimeSeries = SortMask(:,1:iVolume);
        [~, ~, ~, stats] = xASL_stat_ttest(TempTimeSeries,0,0.05,'both',2);
        tValue(iVolume,1) = xASL_stat_MedianNan(stats.tstat(:));
    end
    fprintf('\n');
    
    INDEXn = round(0.5*length(tValue));
    OptimumV = max(tValue(INDEXn:end));
    mintValue = max(find(tValue(INDEXn:end)==OptimumV)+INDEXn-1);
    % max() is added here in coincidental case where there are 2 identical t-values, max() errs on the conservative side
    MinimumtValue = max(tValue(INDEXn:end));
    
    if tValue(end)>(1-x.modules.asl.SpikeRemovalThreshold)*MinimumtValue
        % only exclude frames if the optimal t-value
        % is more than x% higher than including all frames (default
        % x.modules.asl.SpikeRemovalThreshold= 0.01;
        mintValue = length(tValue);
    end
    
    mintValuePlot = zeros(1,length(tValue));
    mintValuePlot(mintValue+1:end)=min(tValue);
    
    % Detect frames for exclusion
    if bSubtraction % if ASL
        exclusion = zeros(1, length(MotionTimeSort)*2);
    else  % if no ASL (e.g. fMRI)
        exclusion = zeros(1, length(MotionTimeSort));
    end
    
    if  mintValue<length(tValue)
        for iFrame=mintValue+1:length(MotionTimeSort)
            % Exclude pair
            ExcludePair = MotionTimeSort(iFrame, 2);
            if  bSubtraction % if ASL
                ExcludeFrames = [ExcludePair*2-1 ExcludePair*2];
            else % if no ASL (e.g. fMRI)
                ExcludeFrames = ExcludePair;
            end
            exclusion(ExcludeFrames) = 1;
        end
    end
    
    if usejava('jvm') % only if JVM loaded
        % Save threshold-free spike detection
        hold on
        subplot(3,1,3);
        hold on
        plot(exclusion,'r');
        hold on
        ylabel('Exclusion matrix');
        axis([1 length(NDV{ii}) 0 ceil(max(NDV{ii}))]); % fix X-axes to be same for subplots
    end
else
    fprintf('ENABLE was skipped...\n');
end

if usejava('jvm') % only if JVM loaded
    jpgfile = fullfile(x.D.MotionDir, ['rp_' x.P.SubjectID '_' x.P.SessionID '_motion.jpg']);
    fprintf('Saving motion plot to %s\n', jpgfile);
    saveas(fig, jpgfile, 'jpg');
    close all;
    clear fig;
end

%only run if ENABLE is on
if bENABLE && bSubtraction && nFrames>10 % if we performed outlier exclusion
    tValue(1:3) = tValue(4); % for nicer plotting
    
    if usejava('jvm') % only if JVM loaded
        fig = figure('Visible','off');
        plot([1:length(tValue)],tValue,'b',[1:length(tValue)],mintValuePlot,'r');
        xlabel('control-label pairs sorted by motion');
        ylabel('mean voxel-wise 1-sample t-test p-value');
        PercExcl    = round((sum(exclusion)/length(exclusion)*100)*10)/10;
        title(['Threshold free motion spike exclusion (red, ' num2str(PercExcl) '%) for ' x.P.SubjectID '_' x.P.SessionID]);
        jpgfile = fullfile( x.D.MotionDir,['rp_' x.P.SubjectID '_' x.P.SessionID '_threshold_free_spike_detection.jpg']);
        fprintf('Saving motion plot to %s\n',jpgfile);
        saveas(fig,jpgfile,'jpg');
        close all;
        clear fig;
    end
    
    % Save 7 images, 3 before & 3 after exclusion
    IndexIs = [1 round(mintValue/3)  round(mintValue/2) mintValue];
    diffIndex = (length(tValue)-mintValue)/3;
    IndexIs(5:7) = [mintValue+diffIndex mintValue+2*diffIndex length(tValue)];
    IndexIs = round(IndexIs);
    Slice2Show = floor(size(IM,3)*0.67); % e.g. slice 11/17
    % pre-allocation for more efficient memory usage
    ExampleIM = zeros(size(SortIM,1), size(SortIM,2), length(IndexIs));
    ExampleIM = single(ExampleIM);
    for iVolume=1:length(IndexIs)
        ExampleIM(:,:,iVolume) = xASL_stat_MeanNan(SortIM(:,:,Slice2Show,1:IndexIs(iVolume)), 4);
    end
    TotalCheck = xASL_vis_TileImages(xASL_im_rotate(ExampleIM,90), 4);
    
    % Find intensities
    SortValues = sort(TotalCheck(isfinite(TotalCheck)));
    MinValue = SortValues(max(1,round(0.001*length(SortValues))));
    MaxValue = SortValues(round(0.999*length(SortValues)));
    
    TotalCheck(TotalCheck<MinValue) = MinValue;
    TotalCheck(TotalCheck>MaxValue) = MaxValue;
    
    jpgfile = fullfile( x.D.MotionDir,['rp_' x.P.SubjectID '_' x.P.SessionID '_PWI_motion_sorted.jpg']);
    fprintf('Saving motion plot to %s\n',jpgfile);
    xASL_vis_Imwrite(TotalCheck, jpgfile);
    
    %% ----------------------------------------------------------------------------------------
    %% 4 Remove spike frames from nifti
    
    fprintf('Remove spike frames from nifti\n');
    
    if sum(exclusion)>0 % only if spikes have been detected
        
        % Load nifti
        TempIm = xASL_io_Nifti2Im(InputPath);
        
        % Remove spikes
        nextFrame = 1;
        for iFrame=1:nFrames
            if ~(exclusion(iFrame))
                NewIm(:,:,:,nextFrame)  = TempIm(:,:,:,iFrame);
                nextFrame = nextFrame+1;
            end
        end
        
        % skip this for fMRI, which is more complicated
        % due to tissue T1 effects (incomplete saturation, so temporal relation between volumes)
        % Save in single precision, since conversion back to INT16 would
        % otherwise loose precious precision
        % PM: can this be simplified to same format as original??
        xASL_io_SaveNifti(x.P.Path_ASL4D,x.P.Path_despiked_ASL4D,NewIm,32,0);
        
        % Do same for *.mat
        LoadParms = load(x.P.Path_ASL4D_mat, '-mat');
        mat = LoadParms.mat;
        
        Incl = [];
        for iExcl=1:length(exclusion)
            if ~exclusion(iExcl)
                Incl(end+1) = iExcl;
            end
        end
        mat = mat(:,:,Incl);
        save(x.P.Path_despiked_ASL4D_mat, 'mat');
    end
    
else
    exclusion = 0;
    PercExcl = 0;
    MinimumtValue = 0;
end

xASL_delete(rInputPath); % delete temporary image

if nFrames>1 % if we performed MoCo
    % Save results for later summarization in analysis module
    save(fullfile(x.D.MotionDir, ['motion_correction_NDV_' x.P.SubjectID '_' x.P.SessionID '.mat']),'NDV','median_NDV','mean_NDV','max_NDV','SD_NDV','MAD_NDV','exclusion','PercExcl','MinimumtValue');
end


end


