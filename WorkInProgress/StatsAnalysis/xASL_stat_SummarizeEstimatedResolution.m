function Summarize_estimated_resolution( x )
%SUMMARIZE_ESTIMATED_RESOLUTION Loads all estimated resolution
% parameters and summarizes them into a Figure


% Skip if these data are not available
if  length(xASL_adm_GetFileList(fullfile(x.D.PopDir,'ResolutionEstimation'), '^.*\.mat$','FPList',[0 Inf]))>=x.dataset.nSubjects

    fprintf('%s\n','Loading files...  '); 
    for iS=1:x.dataset.nSubjects
        xASL_TrackProgress(iS,x.dataset.nSubjects);
        Fname   = fullfile(x.D.PopDir,'ResolutionEstimation',[x.SUBJECTS{iS} '.mat']);
        load(Fname);
        FWHM_all(iS,:)      = optimFWHM_Res_mm;
        Ratio_all(iS,1)     = optimRatio;
    end

    fprintf('\n');
    
    mean_FWHM               = mean(FWHM_all,1);
    SD_FWHM                 = std(FWHM_all,1);

    mean_Ratio              = mean(Ratio_all,1);
    SD_Ratio                = std(Ratio_all,1);

    for ii=1:3
        [N(:,ii) X(:,ii)]       = hist(FWHM_all(:,ii));
    end

    ColorPlot   = {'r' 'b' 'g'};
    fig = figure('Visible','off');

    for ii=1:3
        plot(X(:,ii),N(:,ii),ColorPlot{ii});
        hold on;
    end

    xlabel('FWHM (mm)');
    ylabel('n Subjects');
    title(['Estimated resolution distribution, X (red) Y (blue) Z (green) for GM-WM CBF ratio ' num2str(mean_Ratio) ' +/- ' num2str(SD_Ratio)]);

    OutPutName  = fullfile(x.D.PopDir,'ResolutionEstimation','Distribution_res_est.jpg');

    print(gcf,'-djpeg','-r300', OutPutName   );

end
end
