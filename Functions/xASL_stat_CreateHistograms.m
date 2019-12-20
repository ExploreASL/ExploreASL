function xASL_stat_CreateHistograms( x )
%xASL_stat_CreateHistograms Creates histograms, part of statistical module of ExploreASL
% HJMM Mutsaerts, ExploreASL 2016
%
% Generally, you would permute the datasets (e.g. sessions, cohorts, scanners, etc.) in all possible ways when summarizing and
% comparing the data from set to set using this function.
% Creates single histograms and double, to compare

    %% Initiate
    
    % Repeat for each ROI
    % Start with first dataset, then compare with other dataset if there is any

    plotColorOptions        = {'k'      'g'     'y'         'r'     'c'     'm'         'b'};
    plotColorNames          = {'black'  'green' 'yellow'    'red'   'cyan'  'magenta'   'blue'};

	for iMeas=1:size(x.S.DATASETS{1},2)

            % Define data
            DATA{1}     = x.S.DATASETS{1}(:,iMeas);
            nSubjects   = size(x.S.DATASETS{1},1);

            % Define settings
            % First get range, to scale with
            minMin      = min(DATA{1}(:));
            maxMax      = max(DATA{1}(:));

            rangeRange  = abs(maxMax-minMin);
            xScaler     = rangeRange/5;
            % Get it rounded for power of 10
            xScalerLog10= round(log10(xScaler))-1;
            Rounder     = 10^xScalerLog10;
            xScaler     = round(xScaler/Rounder)*Rounder;



        %% SINGLE/FIRST DATASET
        if length(x.S.DATASETS)==1 % if only single group

            bin_nr      = round(length(DATA{1})/5); % 5 bins, bin size = 1 score, add 1 to lower & upper boundaries to prevent clipping
            min_nr      = round(min(DATA{1}(:)) / xScaler)*xScaler - xScaler; % rounded in 10, 10 extra marge to prevent clipping at edges
            max_nr      = round(max(DATA{1}(:)) / xScaler)*xScaler + xScaler; % same as min but on other side
            bin_size    = (max_nr-min_nr)/bin_nr;
    %         myfilter    = fspecial('gaussian',[bin_nr,1],1*0.00001*bin_nr);

            % Show
            fig = figure('Visible','off');
            [N X]       = hist( DATA{1} (isfinite(DATA{1})) ,[min_nr:(max_nr-min_nr)/bin_nr:max_nr]);
            N           = N./sum(N)./bin_size; % normalization
    %         N           = imfilter(N, myfilter', 'replicate');
            plot(X,100.*N);

            printTitle  = [x.S.Measurements{iMeas} ' ' x.S.output_ID ' histogram ' x.S.NAME ' (n=' num2str(nSubjects) ')'];
            h   = title( printTitle,'interpreter','none'  );
            ylabel('Normalized frequency (a.u.)');
            xlabel([x.S.output_ID ' (' x.S.unit ')']);
            SaveFile   = fullfile( x.S.StatsDir, [printTitle '.jpg']);
            saveas( fig ,SaveFile,'jpg');
            SaveFile   = fullfile( x.S.StatsDir, [printTitle '.eps']);
            saveas( fig ,SaveFile,'epsc');
            
            clear N X bin_nr min_nr max_nr bin_size myfilter printTitle
            close all

        %% COMPARISON WITH SECOND DATASET
        elseif length(x.S.DATASETS)>1 % if another group exists to compare with
%                 DATA{1}             = x.S.DATASETS{1}(:,iMeas);
                DATA{2}             = x.S.DATASETS{2}(:,iMeas);
                nSubjects           = size(x.S.DATASETS{1},1);

                % Define settings
                for ii=1:2
                    bin_nr(ii)      = round(length(DATA{ii})/5); % 5 bins, bin size = 1 score, add 1 to lower & upper boundaries to prevent clipping
                    min_nr(ii)      = round(min(DATA{ii}(:)) / xScaler)*xScaler - xScaler; % rounded in 10, 20 extra marge to prevent clipping at edges
                    max_nr(ii)      = round(max(DATA{ii}(:)) / xScaler)*xScaler + xScaler; % same as min but on other side

                    bin_nr(ii)      = round(length(DATA{ii})/5); % 5 bins, bin size = 1 score, add 1 to lower & upper boundaries to prevent clipping
                    min_nr(ii)      = round(min(DATA{ii}(:)) / xScaler)*xScaler - xScaler; % rounded in 10, 20 extra marge to prevent clipping at edges
                    max_nr(ii)      = round(max(DATA{ii}(:)) / xScaler)*xScaler + xScaler; % same as min but on other side
                end

                bin_nr              = round(mean(bin_nr));
                min_nr              = min(min_nr);
                max_nr              = max(max_nr);

                bin_size            = (max_nr-min_nr)/bin_nr;
        %         myfilter    = fspecial('gaussian',[bin_nr,1],1*0.00001*bin_nr);

                % Show
                fig = figure('Visible','off');
                hold on;

                for ii=1:2
                    [N X]               = hist( DATA{ii} (isfinite(DATA{ii})) ,[min_nr:(max_nr-min_nr)/bin_nr:max_nr]);
                    N                   = N./sum(N)./bin_size; % normalization
            %         N           = imfilter(N, myfilter', 'replicate');
                    plot(X,100.*N,plotColorOptions{ii});
                    clear N X
                end

                printTitle = [x.S.Measurements{iMeas} ' ' x.S.output_ID ' histogram ' x.S.NAME ' (n=' num2str(nSubjects) ')'];

                legend(x.S.COMP{1},S.COMP{2});

                h   = title( printTitle,'interpreter','none'  );
                ylabel('Normalized frequency (a.u.)');
                xlabel([x.S.output_ID ' (' x.S.unit ')']);
                SaveFile   = fullfile( x.S.StatsDir, [printTitle '.jpg']);
                saveas( fig ,SaveFile,'jpg');
                SaveFile   = fullfile( x.S.StatsDir, [printTitle '.eps']);
                saveas( fig ,SaveFile,'epsc');
                
                clear N X bin_nr min_nr max_nr bin_size myfilter printTitle
                close all
        end
    end

end