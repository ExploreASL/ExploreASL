function xASL_stat_PrintBasicStats(x)
%xASL_stat_PrintBasicStats Prints row of basic statistics, part of statistical module of ExploreASL
% HJMM Mutsaerts, ExploreASL 2016
%
% ONE_TWO_SAMPLE_TEST: 1 = one-sample or paired sample
%                      2 = two-sample (unpaired)
% Generally, you would permute the datasets (e.g. sessions, cohorts, scanners, etc.) in all possible ways when summarizing and
% comparing the data from set to set using this function.
 
    %% Initiate
    if ~isfield(x.S,'P_VALUE') || isempty(x.S.P_VALUE)
        x.S.P_VALUE = 0.05;
    end
    
    x.S.SaveFile                 = fullfile(x.S.StatsDir,['stats_p-value_' xASL_num2str(x.S.P_VALUE) '_' x.S.NAME '.tsv']);
 
    % First build cell array, then print to csv-file
    x.S.FID               = fopen(x.S.SaveFile,'wt');
 
    % Print header
    Header1                     = {'Measurement' 'n' 'SW norm H' 'SW norm P' 'mean' 'std' 'CV (%)' 'median' 'MAD' 'CV'};
    for iH=1:length(Header1)
        fprintf(x.S.FID,'%s\t', Header1{iH});
    end
    if length(x.S.DATASETS)>1
        Header2                 = {'n' 'SW norm H' 'SW norm P' 'mean' 'std' 'CV (%)' 'median' 'MAD' 'CV' 't-test H' 't-test P' 'W ranksum H' 'W ranksum P' 'Levene H' 'Levene P' 'B-Forsythe H' 'B-Forsythe P'};
        for iH=1:length(Header2)
            fprintf(x.S.FID,'%s\t', Header2{iH});
        end
    end
    fprintf(x.S.FID,'\n');
 
    fprintf(x.S.FID,'\n%s\n\n',S.NAME);
 
    % Print specific set/condition names
    if length(x.S.DATASETS)>1
        fprintf(x.S.FID,'%s\t\t\t\t\t\t\t\t\t\t%s\n\n',S.COMP{1},S.COMP{2});
 
    end
 
    x.S.ROI     = x.S.Measurements;
 
    % Create single row for each ROI. Start with first dataset, then compare with other dataset if there is any
    for iRoi=1:size(x.S.DATASETS{1},2)
 
        fprintf(x.S.FID,'%s\t', x.S.ROI{iRoi});  % print ROI name
 
        %% SUMMARIZE FIRST DATASET
        DATA1   = x.S.DATASETS{1}(:,iRoi);
 
        fprintf(x.S.FID,'%s\t', xASL_num2str(size(DATA1,1))); % n
 
        % Test for normality
        if  sum(isfinite(DATA1))>2 && size(DATA1,1)<5000 % (restriction valid observations = 3 < n < 5000)
            try 
            [H, P] = xASL_stat_ShapiroWilk(DATA1);         
            catch
                H = NaN;
                P = NaN;
            end
                
        else
            H = NaN;
            P = NaN;
        end
        fprintf(x.S.FID,'%s\t', xASL_num2str(H));
        fprintf(x.S.FID,'%s\t', xASL_num2str(P));
 
        % parametric
        fprintf(x.S.FID,'%s\t', xASL_num2str(xASL_stat_MeanNan(DATA1)));
        fprintf(x.S.FID,'%s\t', xASL_num2str( xASL_stat_StdNan(DATA1,1))); % no degree of freedom correction for n
        fprintf(x.S.FID,'%s\t', xASL_num2str( 100.*(xASL_stat_StdNan(DATA1,1)/xASL_stat_MeanNan(DATA1)) ));
        % non-parametric
        fprintf(x.S.FID,'%s\t', xASL_num2str(xASL_stat_MedianNan(DATA1) ) );
        fprintf(x.S.FID,'%s\t', xASL_num2str(xASL_stat_MadNan(DATA1,1)) ); % Using medians
        fprintf(x.S.FID,'%s\t', xASL_num2str( 100.*(xASL_stat_MadNan(DATA1,1)/xASL_stat_MedianNan(DATA1) ) ));
 
        %% COMPARISON WITH SECOND DATASET
        if length(x.S.DATASETS)>1 % if another group exists to compare with
 
            %% SUMMARIZE 2ND DATASET
            DATA2   = x.S.DATASETS{2}(:,iRoi);
 
            fprintf(x.S.FID,'%s\t', xASL_num2str(size(DATA2,1)));
 
            if  sum(isfinite(DATA1))>2 && size(DATA2,1)<5000
                try
                    [H, P, ~]         = xASL_stat_ShapiroWilk(DATA2);         % (restriction = 3 < n < 5000)
                catch H = NaN;    P = NaN;
                end
            else    H = NaN;    P = NaN;
            end
 
            fprintf(x.S.FID,'%s\t', xASL_num2str(H));
            fprintf(x.S.FID,'%s\t', xASL_num2str(P));
 
            % parametric
            fprintf(x.S.FID,'%s\t', xASL_num2str(xASL_stat_MeanNan(DATA2)));
            fprintf(x.S.FID,'%s\t', xASL_num2str( xASL_stat_StdNan(DATA2,1))); % no degree of freedom correction for n
            fprintf(x.S.FID,'%s\t', xASL_num2str( 100.*(xASL_stat_StdNan(DATA2,1)/xASL_stat_MeanNan(DATA2)) ));
            % non-parametric
            fprintf(x.S.FID,'%s\t', xASL_num2str(xASL_stat_MedianNan(DATA2) ) );
            fprintf(x.S.FID,'%s\t', xASL_num2str(xASL_stat_MadNan(DATA2,1)) ); % Using medians
            fprintf(x.S.FID,'%s\t', xASL_num2str( 100.*(xASL_stat_MadNan(DATA2,1)/xASL_stat_MedianNan(DATA2) ) ));
 
            if ~isfield(x.S,'ONE_TWO_SAMPLE_TEST')
                error('Selection one/paired or two sample unpaired missing!');
            end
 
            % ttest
            if  x.S.ONE_TWO_SAMPLE_TEST==1 % (e.g. for comparing sessions, i.e. within-subject)
                % one sample or paired tests
                % parametric
                if sum(isfinite(DATA1))>2 && size(DATA2,1)>2
                    
                    try 
                        [H, P]                  = xASL_stat_ttest(DATA1,DATA2);
                    catch H = NaN;    P = NaN;
                    end
                else    H = NaN;    P = NaN;
                end
 
                fprintf(x.S.FID,'%s\t', xASL_num2str(H));
                fprintf(x.S.FID,'%s\t', xASL_num2str(P));
 
                % non-parametric
                try  % requires Matlab statistical toolbox
                    [P H]                   = signtest(DATA1,DATA2);
                catch H = NaN;    P = NaN;
                end
                
                fprintf(x.S.FID,'%s\t', xASL_num2str(H));
                fprintf(x.S.FID,'%s\t', xASL_num2str(P));                
 
            elseif  x.S.ONE_TWO_SAMPLE_TEST==2 % (e.g. for comparing groups, i.e. between-subjects)
                % two sample (unpaired) tests
                % parametric
 
                if size(DATA1,1)>2 && size(DATA2,1)>2
                    try 
                        [H, P]                  = xASL_stat_ttest2(DATA1,DATA2);
                    catch H = NaN;    P = NaN;
                    end
                    
                else    H = NaN;    P = NaN;
                end
                
                fprintf(x.S.FID,'%s\t', xASL_num2str(H));
                fprintf(x.S.FID,'%s\t', xASL_num2str(P));
 
 
                % non-parametric
                % two-sample Wilcoxon rank sum test for equal medians
                try % requires Matlab statistical toolbox
                    [P H]               = ranksum(DATA1,DATA2);
                catch H = NaN;    P = NaN;
                end                
                
                fprintf(x.S.FID,'%s\t', xASL_num2str(H));
                fprintf(x.S.FID,'%s\t', xASL_num2str(P));
            end
 
 
            % two-sample Levenetest
            nX                      = length(DATA1);
            nY                      = length(DATA2);
            L_X(1   :nX   ,1)       = DATA1;
            L_X(1   :nX   ,2)       = 1;
            L_X(nX+1:nX+nY,1)       = DATA2;
            L_X(nX+1:nX+nY,2)       = 2;
 
            try 
                [H P]               = xASL_stat_EqualVariancesTest(L_X,0.05,'Levene');
            catch H = NaN;    P = NaN;
            end                      
            
            fprintf(x.S.FID,'%s\t', xASL_num2str(H));
            fprintf(x.S.FID,'%s\t', xASL_num2str(P));
 
            % two-sample Brown-Forsythe test
            try 
                [H P]                   = xASL_stat_EqualVariancesTest(L_X,0.05,'BrownForsythe'); 
            catch H = NaN;    P = NaN;
            end                      
            
            fprintf(x.S.FID,'%s\t', xASL_num2str(H));
            fprintf(x.S.FID,'%s\t', xASL_num2str(P));
        end
        fprintf(x.S.FID,'\n');
    end
 
    fclose(x.S.FID);
 
 
end
