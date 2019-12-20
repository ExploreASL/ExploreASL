%% xASL_create_wsCV_map_reproducibility

Bmask       = xASL_io_Nifti2Im('C:\ExploreASL\Maps\WBmaskASL.nii');
SaveDir     = 'C:\Backup\ASL\xASL_paper_repeat_analysis\AnalysisPrints';

% Load data

Ana{1}      = {'analysis-1_sub-HC_Q-0_OS-Lin' 'analysis-2_sub-HC_Q-0_OS-Lin'};
Ana{2}      = {'analysis-1_sub-HC_Q-1_OS-Lin' 'analysis-2_sub-HC_Q-1_OS-Lin'};
Ana{3}      = {'analysis-1_sub-HC_Q-0_OS-Win' 'analysis-2_sub-HC_Q-0_OS-Win'};
Ana{4}      = {'analysis-1_sub-HC_Q-1_OS-Win' 'analysis-2_sub-HC_Q-1_OS-Win'};
Ana{5}      = {'analysis-1_sub-HC_Q-0_OS-Win' 'analysis-1_sub-HC_Q-0_OS-Lin'};
Ana{6}      = {'analysis-1_sub-HC_Q-1_OS-Win' 'analysis-1_sub-HC_Q-1_OS-Lin'};

for iA=1:length(Ana)
    ROOT{1}     = fullfile('C:\Backup\ASL\xASL_paper_repeat_analysis',Ana{iA}{1},'dartel\Stats');
    ROOT{2}     = fullfile('C:\Backup\ASL\xASL_paper_repeat_analysis',Ana{iA}{2},'dartel\Stats');

    clear RegExp
    RegExp{1}   = '^spatialCoV_qCBF_untreated_TotalGM_.*PVC2\.csv$';
    RegExp{2}   = '^mean_qCBF_untreated_TotalGM_.*PVC2\.csv$';
    RegExp{3}   = '^mean_qCBF_untreated_DeepWM_.*PVC2\.csv$';

    for iR=1:2
        for iT=1:3
            FilePath{iT}        = xASL_adm_GetFileList(ROOT{iR},RegExp{iT});
            [numericData, textData, rawData] = xlsread(FilePath{iT}{1});

            ReproData{iT}(:,iR)      = cell2mat(rawData(3:end,5));
        end
    end

    %% Plot PVC reproducibility values

    LabelNames     = {'Spatial CoV (norm PVC, ratio)' 'GM CBF (PVC, mL/100g/min)' 'Deep WM CBF (PVC, mL/100g/min)'};

    for iT=1:3
        % Compute values
        Mean2       = mean(ReproData{iT},2);
        Diff2       = ReproData{iT}(:,1) - ReproData{iT}(:,2);
        MeanAll     = mean(Mean2);
        MeanDiff    = mean(Diff2);
        SDdiff      = std(Diff2);
        wsCV        = (SDdiff/MeanAll)*100;
        % Create plot
        close all
        figure(iT);
        plot(Mean2,Diff2,'k.'); % data
        hold on
        plot([min(Mean2) max(Mean2)],repmat(MeanDiff,2),'k-'); % data
        hold on
        plot([min(Mean2) max(Mean2)],repmat(MeanDiff-1.96*SDdiff,2),'k--'); % data
        hold on
        plot([min(Mean2) max(Mean2)],repmat(MeanDiff+1.96*SDdiff,2),'k--'); % data

        xlabel(['Mean ' LabelNames{iT}]);
        ylabel(['Diff ' LabelNames{iT}]);
        PrintName   = fullfile(SaveDir, [num2str(iA) '_' Ana{iA}{1}(19:end) '_' LabelNames{iT}(1:13)]);
        if  exist(PrintName,'file')
            delete(PrintName);
        end
        print(PrintName,'-djpeg');
    end

    %% Show wsCV images

    % Load data
    ROOT{1}     = fullfile('C:\Backup\ASL\xASL_paper_repeat_analysis',Ana{iA}{1},'dartel');
    ROOT{2}     = fullfile('C:\Backup\ASL\xASL_paper_repeat_analysis',Ana{iA}{2},'dartel');

    RegExp      = '^qCBF_untreated_.*\.(nii|nii\.gz)$';

    for iR=1:2
        IMlist        = xASL_adm_GetFileList(ROOT{iR},RegExp);
        for iM=1:length(IMlist)
            xASL_TrackProgress(iM,length(IMlist));
            IM(:,:,:,iM,iR)    = xASL_io_Nifti2Im(IMlist{iM});
        end
    end

    % Compute maps
    Mean2       = xASL_stat_MeanNan(IM,5);
    Diff2       = IM(:,:,:,:,1) - IM(:,:,:,:,2);
    MeanAll     = xASL_stat_MeanNan(Mean2,4);
    MeanDiff    = xASL_stat_MeanNan(Diff2,4);
    SDdiff      = xASL_stat_StdNan(Diff2,[],4);
    wsCV        = (SDdiff./MeanAll)*100;
    wsCV        = wsCV.*logical(Bmask);

    wsCV_im     = TransformDataViewDimension(wsCV);
    jet_256     = jet(256); jet_256(1,:) = 0;
    close all

    if iA<5
        MaxInt  = 0.25;
    elseif iA==5
        MaxInt  = 10;
    else
        MaxInt  = 2.5;
    end

    figure(1);imshow(wsCV_im,[0 MaxInt],'colormap',jet_256,'border','tight')
    colorbar;
    PrintName   = fullfile(SaveDir, [num2str(iA) '_' Ana{iA}{1}(19:end) '_wsCV_im.jpg']);
    print(PrintName,'-djpeg');
end
