function xASL_wrp_PermuteWithinSet( x )
%xASL_wrp_PermuteWithinSet Part of statistical module of ExploreASL
% HJMM Mutsaerts, ExploreASL 2016
% INPUT x.S structure should contain fields:
% populations_mat
% pop_name
% ROI
% nROI

    warning('off','images:initSize:adjustingMag'); % warning about scaling if an image doesnt fit screen, disable

    if ~isfield(x.S, 'ONE_TWO_SAMPLE_TEST')
        error('Selection one/paired or two sample unpaired missing!');
    end
 
    x.S.nSets              = length(x.S.DATASETS_RESTR);
    done{1}                = 'NaN'; % start keeping track
    OriginalStatsDir       = x.S.StatsDir;
    
    %% 1 sample t-test is often not interesting, only run it in case of co-variates, for regression
    HasCovariates           = max(x.S.Sets1_2Sample==3);

    if  HasCovariates
        for iSet=1:x.S.nSets
    %             x.S.StatsDir      = fullfile( OriginalStatsDir, ['Set_' num2str(iSet)]);
    %             xASL_adm_CreateDir(x.S.StatsDir);
            x.S.StatsDir      = fullfile( OriginalStatsDir);
            x.S.NAME          = [x.S.SetsName{x.S.iCurrentSet} x.S.SetsOptions{x.S.iCurrentSet}{iSet}];
            x.S.DATASETS{1}   = x.S.DATASETS_RESTR{iSet}; % specify data
            x.S.CoVar{1}      = x.S.CoVar_RESTR{iSet}; % specify data
            x.S.COMP{1}       = x.S.SetsOptions{x.S.iCurrentSet}{iSet};
            x.S.iSet{1}       = iSet;

            SaveData(x);
%             x.S.function2call( x ); % function 2 call (e.g. ROI_PRINT_STATS or hist)    
            x.S               = rmfield(x.S, {'iSet' 'DATASETS' 'COMP' 'NAME' 'CoVar'});        
        end
    end
    
    x.S.StatsDir          = OriginalStatsDir;    

%     if x.S.KISS~=3 % 3 requests to skip the statistics, but only permute sets themselves to show descriptives
    
        %% Do ANOVA first if multiple sets
        if x.S.nSets>2
            x.S.NAME          = x.S.SetsName{x.S.iCurrentSet};
            for iSet=1:x.S.nSets
                x.S.NAME          = [x.S.NAME ' ~ ' x.S.SetsOptions{x.S.iCurrentSet}{iSet}];
                x.S.COMP{iSet}    = x.S.SetsOptions{x.S.iCurrentSet}{iSet};
                x.S.iSet{iSet}    = iSet;
            end

            x.S.DATASETS          = x.S.DATASETS_RESTR; % simply use all data in this case
            x.S.CoVar             = x.S.CoVar_RESTR; % simply use all data in this case
            
            if isfield(x.S,'PerSetName')
                x.S.NAME      = [x.S.NAME ' ' x.S.PerSetName];
            end    

            SaveData(x);           
            % x.S.function2call( x ); % function 2 call (e.g. ROI_PRINT_STATS or hist)    
            x.S               = rmfield(x.S, {'iSet' 'DATASETS' 'COMP' 'NAME' 'CoVar'});
        end


        %% ITERATE SESSIONS (e.g. t-tests)
        for iComp=1:x.S.nSets % Population to be compared against
            x.S.iSet{1}       = iComp;
           for iSumm=1:x.S.nSets % Population to be summarized
                x.S.iSet{2}   = iSumm;

                % Keeping track of which comparisons have been made, to avoid illustrating double
                % Session summary will be printed in same row as comparison with another session
                % Since we allow sessions to be compared against themselves, 
                % there will never be a session or comparison skipped in this way.

                performed=0;

                for ii=1:length(done) % for each performed combination
                    if ~isempty(strfind( done{ii},[num2str(iComp) '_' num2str(iSumm)] )) || ~isempty(strfind( done{ii},[num2str(iSumm) '_' num2str(iComp)] ))
                        performed=1; % if any of the permutations of the to be summarized and compared against sessions exists, skip iteration
                    end
                end

                % if it hasn't been performed yet, and there are two different datasets, let's do it!
                %%% -> if there are no "sets/conditions", then "all scans" will correctly summarize the data
                %%% -> conditions should not be summarized against theirselves, this causes redundant printing
                if ~performed && iComp~=iSumm


                    x.S.NAME          = [x.S.SetsName{x.S.iCurrentSet} ' ' x.S.SetsOptions{x.S.iCurrentSet}{iSumm} ' vs. ' x.S.SetsOptions{x.S.iCurrentSet}{iComp}];

                    if isfield(x.S,'PerSetName')
                        x.S.NAME      = [x.S.NAME ' ' x.S.PerSetName];
                    end

                    x.S.DATASETS{1}   = x.S.DATASETS_RESTR{iComp}; % specify data
                    x.S.DATASETS{2}   = x.S.DATASETS_RESTR{iSumm}; % specify data

                    x.S.CoVar{1}      = x.S.CoVar_RESTR{iComp}; % specify covariants
                    x.S.CoVar{2}      = x.S.CoVar_RESTR{iSumm}; % specify covariants
                    
                    x.S.COMP{1}       = x.S.SetsOptions{x.S.iCurrentSet}{iComp};
                    x.S.COMP{2}       = x.S.SetsOptions{x.S.iCurrentSet}{iSumm};

                    x.S.iSet{1}       = iComp;
                    x.S.iSet{2}       = iSumm;                

 
                    SaveData(x);
%                     x.S.function2call( x ); % function 2 call (e.g. ROI_PRINT_STATS or hist)

                    % Keeping track: note that this one has been performed
                    done{length(done)+1}   = [num2str(iComp) '_' num2str(iSumm)];

                    % Remove residual fields
                    x.S               = rmfield(x.S, {'DATASETS' 'COMP' 'NAME' 'CoVar'});
                end
            end
        end
%     end
end
 



%% ----------------------------------------------------------------------------------------
%% Save the permuted data

function SaveData(x)
%SaveData Store the data required for the GLM

    % Save data
    if  isfield(x.S,'DAT')
        x.S                         = rmfield(x.S,'DAT');    
    end
    if  isfield(x.S,'DATASETS_RESTR')
        x.S                         = rmfield(x.S,'DATASETS_RESTR');
    end    
    if  isfield(x.S,'CoVar_RESTR')
        x.S                         = rmfield(x.S,'CoVar_RESTR');
    end        
    
    xMatFile                    = fullfile(x.S.StatsDir,['xASL_' x.S.NAME '.mat']);
    fprintf('%s\n',['Saving ' xMatFile]);
    save(xMatFile,'x');

end



