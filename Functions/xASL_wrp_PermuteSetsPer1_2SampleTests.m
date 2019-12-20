function xASL_wrp_PermuteSetsPer1_2SampleTests( x )
%xASL_wrp_PermuteSetsPer1_2SampleTests 
% INPUT
% x = from ExploreASL
% masks   = physically loaded masks for each subject
% ASL     = the data to be analyzed. Could be ASL, or e.g. SD or SNR masks
%
%
% By HJMM Mutsaerts, ExploreASL 2016
%
% The idea behind this permutation function is that sometimes you want
% to treat 1 or 2 sample datasets differently.
% E.g. when doing a spaghetti plot, you want the 2 sample datasets for different colors,
% and the 1 sample datasets on the x-axis

% Requires fields:
% x.S.SetsID
% x.S.SetsName
% x.S.SetsOptions
% x.S.Sets1_2Sample

%% Administration

    if      size(x.S.SetsID,2)~=length(x.S.SetsName)
            error('size(x.S.SetsID,2)~=length(x.S.SetsName)');
    elseif  size(x.S.SetsID,2)~=length(x.S.SetsOptions)
            error('size(x.S.SetsID,2)~=length(x.S.SetsOptions)');
    elseif  size(x.S.SetsID,2)~=length(x.S.Sets1_2Sample)
            error('size(x.S.SetsID,2)~=length(x.S.Sets1_2Sample)');
    end

%% Run    
    for iSet2=1:length(x.S.SetsName)         % Loop sets to get a 2 sample set
        if x.S.Sets1_2Sample(iSet2)==2

            x.S.NAME_SAMPLE_2                 = x.S.SetsName{iSet2};
            x.S.OPTIONS_SAMPLE_2              = x.S.SetsOptions{iSet2};

            for iSet1=1:length(x.S.SetsName) % Loop sets to get a 1 sample set
                if x.S.Sets1_2Sample(iSet1)==1

                    x.S.NAME_SAMPLE_1         = x.S.SetsName{iSet1};
                    x.S.OPTIONS_SAMPLE_1      = x.S.SetsOptions{iSet1};

                    % Restructure for 2 sample set first
                    code_column             = x.S.SetsID(:,iSet2);     % set  column
                    data_columns            = x.S.DAT;                % data columns
                    x.S.DATASETS              = xASL_stat_RestructureSets( data_columns, code_column );
                    clear data_columns

                    % Do the same with the code of the 1 sample dataset (setIDs)
                    data_columns            = x.S.SetsID(:,iSet1);
                    code_columns_new        = xASL_stat_RestructureSets( data_columns, code_column );
                    clear code_column data_columns

                    % Restructure the DATASETS_RESTR into the 1 sample dataset subgroups
                    for iSub=1:length(x.S.DATASETS)
                        data_columns            = x.S.DATASETS{iSub};
                        code_column             = code_columns_new{iSub};
                        x.S.DATASETS_RESTR{iSub}  = xASL_stat_RestructureSets( data_columns, code_column );
                    end

                    x.S.function2call( x ); % function 2 call (e.g. spaghetti plot)
                end
            end
        end
    end



 
 end