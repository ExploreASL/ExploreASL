function [ StatCalc AUC  ] = GetRocCurve( x, GoldPos, FigureOut, FilterFile )
%GetRocCurve

    FilterFile          = xASL_io_ReadNifti( FilterFile );
    FilterFile          = FilterFile.dat(:,:,:);

    meanGMmask          = fullfile( x.D.PopDir, 'DARTEL_T1_template.nii');
    meanGMmask          = nifti ( meanGMmask );
    meanGMmask          = single( meanGMmask.dat(:,:,:) );
    meanGMmask          = meanGMmask./max(meanGMmask(:));
    meanGMmask          = meanGMmask>0.5;

    tMin            = 0;   % minimal t-value (statistical threshold)
    tMax            = max(FilterFile(:));  % maximal t-value (statistical threshold)
    ROCRes          = 100;   % resolution, or nSteps
    StepSize        = (tMax - tMin) ./ ROCRes;

    StatCalc(:,1)   = [tMin:StepSize:tMax]';

    GoldPos         = GoldPos .* meanGMmask(:,:,:,1);

    for iStat       = 1:size(StatCalc(:,1))
        H           = FilterFile>StatCalc(iStat,1);
        [StatCalc(iStat,2) StatCalc(iStat,3)]   = CalcROC( GoldPos, H);
    end

    ind1    = find(StatCalc(:,2)==1);
    if  isempty(find(StatCalc(ind1,3)==1)) % if there is no closing end at [1 1]
        NewStatCalc                             = zeros(size(StatCalc,1),size(StatCalc,2) );
        NewStatCalc(2:size(StatCalc,1)+1,:)     = StatCalc;
        NewStatCalc(1,:)                        = [0 1 1];
        StatCalc                                = NewStatCalc;
    end

    ind1    = find(StatCalc(:,2)==0);
    if  isempty(find(StatCalc(ind1,3)==0)) % if there is no closing end at [0 0]
        NewStatCalc                             = zeros(size(StatCalc,1),size(StatCalc,2) );
        NewStatCalc(2:size(StatCalc,1)+1,:)     = StatCalc;
        NewStatCalc(1,:)                        = [0 1 1];
        StatCalc                                = NewStatCalc;
    end

    if      exist('FigureOut','var')
            fig = figure('Visible','off');
    else    figure(1);
    end

    plot(StatCalc(:,3),StatCalc(:,2)); % plot ROC
    hold on
    plot([0:0.01:1],[0:0.01:1],'k--'); % plot line of equality
    xlabel('False positive rate (1-specificity)');
    ylabel('True positive rate');
    AUC     = trapz(StatCalc(:,3),StatCalc(:,2));
    AUC     = abs(round(AUC*1000)/1000); % round to 3 floating numbers
    text(0.7,0.2,['AUC=' num2str(AUC)]);

    if  exist('FigureOut','var')

        [Path2 file ext]             = fileparts( FigureOut );
        xASL_adm_CreateDir( Path2 );

        print(gcf,'-djpeg','-r200', FigureOut );

    end

end
