function xASL_wrp_PermuteOverSets( x )
%xASL_wrp_PermuteOverSets Permute data across sets (e.g. groups, scanners, cohorts, conditions etc)
% HJMM Mutsaerts, ExploreASL 2016
%
% x.S.function2call   = function to run
% x.S.SetsName        = cell array with different (data-)sets
% x.S.SetsOptions     = names within sets (e.g. before/after medication, scanner-name, session 1`, diseased/controls etc.)

% Lay-out:
% 1 run function for total data
% 2 for all sets:
% A) without any other splitting, permute across this individual set (e.g. session 1 vs. 2, 1 vs. 3, 2 vs. 3
% B) split other sets and permute again across this individual set (e.g. for the diseased only, compare session 1 vs. 2, 1 vs. 3 & 2 vs. 3)

% Within total data & within splitted data, run regressions with any
% available continuous variables.

% Lay-out subfunctions:
% 1) single dataset (e.g. total data or single splitted datasets) will run a
% 1-sample t-test
% 2) two datasets will start paired or two-sample t-tests
% 3) multiple datasets will start 1-way ANOVA
% 4) single dataset with continuous variable will start multiple regression


%% THIS FUNCTION SPLITS DATASETS, other function permutes across them

%% 0    Administration
    % Validity check
    % Everywhere, the same number of sets should be specified

    if isfield(x.S,'SetsName')
        if      length(x.S.SetsName)~=length(x.S.SetsOptions)
                error('Nr of setsnames ~= nr of setsoptions!');
        elseif  length(x.S.SetsName)~=size(x.S.SetsID,2)
                error('Nr of setsnames ~= size(x.S.SetsID,2)!');
        elseif  length(x.S.SetsName)~=length(x.S.Sets1_2Sample)
                error('Nr of setsnames ~= length(x.S.Sets1_2Sample)!');
        end
    end

    if  isfield(x.S,'DATASETS_RESTR')
        x.S       = rmfield(x.S,'DATASETS_RESTR');
    end
    if  isfield(x.S,'CoVar_RESTR')
        x.S       = rmfield(x.S,'CoVar_RESTR');
    end    
    
 
%% 1    Run function for total data

%     if  x.S.KISS==0
        x.S.NAME                      = 'AllScans';

        x.S.StatsDir                  = fullfile( x.S.oriDIR, x.S.NAME );
        xASL_adm_CreateDir(x.S.StatsDir);
        x.S.DATASETS                  = {x.S.DAT};              % specify data, use all data as DATASET
        x.S.CoVar                     = x.S.SetsID ;
        
        SaveData(x);
%         x.S.function2call( x );

        x.S                           = rmfield(x.S,{'DATASETS' 'NAME' 'StatsDir' 'CoVar'});
%     end
   
    % Only perform permutation if there are sets/conditions
    if  isfield(x.S,'SetsName')
        fprintf('%s\n','Performing stats...  ');
        for iSet=1:length(x.S.SetsName)
            xASL_TrackProgress(iSet,length(x.S.SetsName));
            if x.S.Sets1_2Sample(iSet)~=3 && length(x.S.SetsOptions{iSet})>1 && length(x.S.SetsOptions{iSet})<10 
                % if set is not continuous, and has >1 options, but not too many
                
                if  size(x.S.DAT,1)==size(x.S.SetsID,1)  
                    % this is a requirement for the restructuring
                    % e.g. skip between-session comparisons if analyzing
                    % volume or TT that exists only once per subject

                %% 2    Permute across individual set

                    x.S.ONE_TWO_SAMPLE_TEST   = x.S.Sets1_2Sample(iSet);
                    x.S.iCurrentSet           = iSet;
                    
                    SufficientGroupSize     = CheckSufficientSize( x.S.SetsID(:,iSet), 16 ); % Check group size
                    % Here we only continue when we have n=16 for each group
                    IsContinue              = SufficientGroupSize;
                    
%                     if  x.S.KISS==3 && length(unique(x.S.SetsID(:,iSet)))==1
%                         % if we chose to keep it simple, to only do descriptives (x.S.KISS==3) and this set only
%                         % as 1 option (nothing to vary across), we will
%                         % skip this
%                         
%                         IsContinue  = 0;
%                     end
                    
                    
                    if  IsContinue
                        % for each group we split to, we need at least 2
                        % SubjectSessions (e.g. t-test requires at least 2)
                        
                        x.S.StatsDir              = fullfile( x.S.oriDIR, ['permute_' x.S.SetsName{iSet}]);
                        xASL_adm_CreateDir(x.S.StatsDir);
                        
                        cd(x.S.StatsDir);
                        
                        % Restructure sets
                        x.S.CoVar_RESTR           = xASL_stat_RestructureSets( x.S.SetsID, x.S.SetsID(:,iSet) );
                        x.S.DATASETS_RESTR        = xASL_stat_RestructureSets( x.S.DAT, x.S.SetsID(:,iSet) );
                        
                        
                        xASL_wrp_PermuteWithinSet( x ); % permutes within single set (such as sessions)
                        x.S                       = rmfield(x.S,{'DATASETS_RESTR' 'StatsDir' 'CoVar_RESTR'});

                %% 3    Split other ordinal (non-continuous, but representing conditions such as cohort or session),
                %       sets & permute again across this set
                        % Define other sets

%                         if  x.S.KISS~=3 % 3 requests to skip the statistics, but only permute sets themselves to show descriptives
                        
                            if  length(x.S.SetsName)>1 % if there are multiple sets
                                next = 1;
                                for ii=1:length(x.S.SetsName)
                                    if ii~=iSet && x.S.Sets1_2Sample(ii)~=3 && length(x.S.SetsOptions{ii})>1 && length(x.S.SetsOptions{ii})<10
                                        % if set is not the same as CurrentSet, not continuous, & has >1, but not too many, options
                                    otherSets(next)     = ii;
                                    next                = next+1;
                                    end
                                end

                                if  exist('otherSets','var') % if there are other non-continuous sets

                                    for othersetn=1:length(otherSets)
                                        iOtherSet           = otherSets(othersetn);

                                        if ~( strcmp(x.S.OutputID,'volume') && strcmp(x.S.SetsName(iOtherSet),'session') ) % skip session comparison if analyzing volume
                                            % This should be improved by a more general
                                            % checking whether the parameter under
                                            % comparison varies (e.g. volume doesn't
                                            % vary with session)

                                            % restructure data in current otherSet
                                            CodeColumn            = x.S.SetsID(:,iOtherSet);    % set  column

                                            [SufficientGroupSize] = CheckSufficientSize( CodeColumn, 4 ); % Check group size

                                            if  SufficientGroupSize % if sufficiently large groups, split them & run stats on subgroups

                                %                 data_columns        = x.S.DAT;                    
                                                DATASETS            = xASL_stat_RestructureSets( x.S.DAT,    CodeColumn );
                                                CoVar_RESTR         = xASL_stat_RestructureSets( x.S.SetsID, CodeColumn );
                                                DataColumns         = x.S.SetsID(:,iSet);         % do the same for current Set ID
                                                CodeColumnNew       = xASL_stat_RestructureSets( DataColumns, CodeColumn );

                                                x.S.nSets2            = length(x.S.SetsOptions{iOtherSet}); % number of conditions in "OtherSet"

                                                % Iterate within current otherSet to compare within current set per current otherSet

                                                for iSet2=1:x.S.nSets2 % iterate otherSets to compare currentSet per otherSet
                                                    if  x.S.Sets1_2Sample(iOtherSet)~=3

                                                        CodeColumn         = CodeColumnNew{iSet2};       % session column
                                                        DataColumns        = DATASETS{iSet2};              % data columns
                                                        CoVarColumns       = CoVar_RESTR{iSet2};

                                                        if  length(unique(CodeColumn))==length(x.S.SetsOptions{iSet})
                                                            % If not all setsoptions are left when splitting
                                                            % up the other sets, this part will be omitted

                                                            [SufficientGroupSize] = CheckSufficientSize( CodeColumn, 4 ); % Check group size                                                        

                                                            if  SufficientGroupSize

                                                                x.S.DATASETS_RESTR    = xASL_stat_RestructureSets( DataColumns , CodeColumn );
                                                                x.S.CoVar_RESTR       = xASL_stat_RestructureSets( CoVarColumns, CodeColumn );                                

                                                                x.S.PerSetName        = ['within ' x.S.SetsName{iOtherSet} ' ' x.S.SetsOptions{iOtherSet}{iSet2}];
                                                                x.S.StatsDir          = fullfile( x.S.oriDIR, ['permute_' x.S.SetsName{iSet} ' ' x.S.PerSetName]);
                                                                xASL_adm_CreateDir(x.S.StatsDir);

                                                                xASL_wrp_PermuteWithinSet( x );

                                                                % Remove residual fields
                                                                x.S                   = rmfield(x.S, {'DATASETS_RESTR' 'StatsDir' 'PerSetName' 'CoVar_RESTR'});
                                                            end
                                                        else
                                                            x.S.PerSetName        = ['within ' x.S.SetsName{iOtherSet} ' ' x.S.SetsOptions{iOtherSet}{iSet2} '_too_small_sample_size'];
                                                            x.S.StatsDir          = fullfile( x.S.oriDIR, ['permute_' x.S.SetsName{iSet} ' ' x.S.PerSetName]);
                                                            xASL_adm_CreateDir(x.S.StatsDir);                                                
                                                            x.S                   = rmfield(x.S, {'StatsDir' 'PerSetName'});
                                                        end

                                                    end
                                                end

                                            % Remove residual fields
                                            x.S                   = rmfield(x.S, {'nSets2'});

                                            end
                                        end
                                    end
                                end
                            end % Split other
%                         end
                    end
                end
            end
        end
        fprintf('\n');
    end
end % permute over sets

function [SufficientGroupSize] = CheckSufficientSize( CodeColumn, MinimalNdatasets )
    % CheckSufficientSize Checks whether size of groups is sufficiently large
    % to split them for subgroup stats
    UniqueCode          = unique( CodeColumn(isfinite(CodeColumn) ));
    SufficientGroupSize = 1;
    for iU=1:length(UniqueCode)
        if  sum(CodeColumn==UniqueCode(iU))<MinimalNdatasets
            % There are code options that have too few subjects or
            % sessions, so regrouping will go wrong
            SufficientGroupSize =0;
        end
    end
    
    if  max(CodeColumn)~=length(UniqueCode)
        SufficientGroupSize =0;
        % These codes are not indices/options (e.g. 1 2 3 4 5)
        % but rather continuous values
    end
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
