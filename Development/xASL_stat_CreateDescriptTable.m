% xASL_stat_CreateDescriptTable(x);
%xASL_stat_CreateDescriptTable Summary of this function goes here
%   Detailed explanation goes here


x.S.StatsDir % folder where we will store the Table

%% Here, there is space to create customary sets/parameters
%  Create normative WM volumetrics
WMset   = 0;
WMHset  = 0;
for iSet=1:length(x.S.SetsName)
    if      strcmp(x.S.SetsName{iSet},'WM_vol')
            WMset           = iSet;
    elseif  strcmp(x.S.SetsName{iSet},'WMH_vol')
            WMHset          = iSet;
    end
end
if  (WMset*WMHset)~=0 % we have both sets, so we can calculate NAWM_WMRatio & WMH_WMRatio
    x.S.SetsName{end+1}     = 'NAWM_WMRatio';
    x.S.SetsOptions{end+1}  = 'NAWM_WMRatio';
    x.S.Sets1_2Sample(end+1)= 3;
    
    WM_vol                  = x.S.SetsID(:,WMset);
    WMH_vol                 = x.S.SetsID(:,WMHset);
    NAWM_vol                = WM_vol-WMH_vol;
    
    x.S.SetsID(:,end+1)     = NAWM_vol./WM_vol;
    x.S.SetsName{end+1}     = 'WMH_WMRatio';
    x.S.SetsOptions{end+1}  = 'NAWM_WMRatio';
    x.S.SetsID(:,end+1)     = WMH_vol./WM_vol;
    x.S.Sets1_2Sample(end+1)= 3;
end
    
    

%% First, we create a list of data we want to include
% Here, define the sets that we want to include, in the order of appearance in the Table.
% A set will only be included if it has more than 1 unique values
IncludedSets    = {'Age' 'GM_ICVRatio' 'NAWM_WMRatio' 'WMH_WMRatio' 'WMH_count' 'MeanMotion' 'CBF' 'spatial_CoV'};
StatsSumm       = {'Par' 'Par'         'Par'          'Par'         'Par'       'Par'        'Par' 'Par'};
TableSets       = [];

for iInc=1:length(IncludedSets)
    for iSet=1:length(x.S.SetsName)
        IsSameSet   = strcmp(lower(IncludedSets{iInc}),lower(x.S.SetsName{iSet}));
        IsEnoughN   = length(unique(x.S.SetsID(:,iSet)))>1;
        if  IsSameSet && IsEnoughN
            TableSets(end+1)   = iSet;
        end
    end
end

%% Set we want to stratify over
%  If this is empty, it will be ignored
StratifSet      = 'Cohort';
for iSet=1:length(x.S.SetsName)
    if  strcmp(x.S.SetsName{iSet},'Cohort')
        iStratSet     = iSet;
    end
end

%% Loop across each set
clear CurrName MeanSet SDset TableOut StratCurrSetData
for iT=1:length(TableSets)
    % Determine name
    iSet            = TableSets(iT);
    CurrName{iT}    = x.S.SetsName{iSet};
    % Determine unit
    if     ~isempty(strfind(lower(CurrName{iT}),'vol'))
            CurrUnit{iT}  = 'mL';
    elseif ~isempty(strfind(lower(CurrName{iT}),'ratio'))
            CurrUnit{iT}  = '%';
    elseif ~isempty(strfind(lower(CurrName{iT}),'age'))
            CurrUnit{iT}  = 'years';
    elseif ~isempty(strfind(lower(CurrName{iT}),'count'))
            CurrUnit{iT}  = '';
    elseif ~isempty(strfind(lower(CurrName{iT}),'motion'))
            CurrUnit{iT}  = 'mm';
    elseif ~isempty(strfind(lower(CurrName{iT}),'time'))
            CurrUnit{iT}  = 'hh:mm';
    else
            CurrUnit{iT}  = [];
    end
    
    if ~isempty(CurrUnit{iT}) % combine parameter name & unit
        CurrName{iT}  = [CurrName{iT} ' (' CurrUnit{iT} ')'];
    end
    
    switch StatsSumm{iT}
        
        case 'Par' % parametric description, e.g. age=continuous
            CurrSetData     =   x.S.SetsID(:,iSet);
            % Check here if we need to use values from options
            if  x.S.Sets1_2Sample(CurrSet)==2 && ~isempty(str2num(x.S.SetsOptions{iSet}{1})) && ~isempty(str2num(x.S.SetsOptions{iSet}{2}))
                UniqueV     = unique(CurrSetData);
                for iU=1:length(UniqueV)
                    CurrSetData(CurrSetData==UniqueV(iU)) = str2num(x.S.SetsOptions{iSet}{iU});
                end
            end

            GroupN(1)       = numel(CurrSetData);
            GroupName{1}    = 'All';
            MeanSet(iT,1)   = mean(CurrSetData); % could be replaced by median
            SDset(iT,1)     = std(CurrSetData);  % could be replaced by mad
            
            % Now do the same for the stratified Table
            if  exist('iStratSet','var')
                StratData   = x.S.SetsID(:,iStratSet);
                UniqStrat   = unique(StratData);
                for iU=1:length(UniqStrat)
                    StratIncl               = StratData==UniqStrat(iU);
                    StratCurrSetData{iU}    = CurrSetData(StratIncl);
                    GroupName{1+iU}         = x.S.SetsOptions{iStratSet}{iU};
                    GroupN(1+iU)            = numel(StratCurrSetData{iU});
                    MeanSet(iT,1+iU)        = mean(StratCurrSetData{iU});
                    SDset(iT,1+iU)          = std(StratCurrSetData{iU});
                end
                % Add here t-test or ANOVA
                [h,p(iT),ci,stats]          = ttest2(StratCurrSetData{1},StratCurrSetData{2});
                % [p(iT),stats]               = ranksum(StratCurrSetData{1},StratCurrSetData{2});
            end
            
            if  strcmp(CurrUnit{iT},'%') % for percentages, multiply by 100
                MeanSet(iT,:)   = MeanSet(iT,:).*100;
                SDset(iT,:)     = SDset(iT,:).*100;
            end                    
                    
        case 'Cat' % categorical description, e.g. sex: M/F
            % if there are only two possibilities (e.g. Male/Female)
            % we put them in the same Table cell. Otherwise, we make multiple rows from them
        otherwise
    end
end 

%% Create the Table
if  exist('iStratSet','var')
    for iT=1:length(CurrName)
        TableOut{1,1}     = 'Cohort';
        TableOut{1,2}     = [GroupName{2}  ', n=' num2str(GroupN(2)) ' (' num2str(round(GroupN(2)/GroupN(1)*100,1)) '%)'];
        TableOut{1,3}     = [GroupName{3}  ', n=' num2str(GroupN(3)) ' (' num2str(round(GroupN(3)/GroupN(1)*100,1)) '%)'];
        TableOut{1,4}     = 'T-test (p-value)';
        % TableOut{1,4}     = [GroupName{1}  ', n=' num2str(GroupN(1))];
        
        TableOut{1+iT,1}  = CurrName{iT};
        TableOut{1+iT,2}  = [num2str(round(MeanSet(iT,2),1)) ' ? ' num2str(round(SDset(iT,2),1))];
        TableOut{1+iT,3}  = [num2str(round(MeanSet(iT,3),1)) ' ? ' num2str(round(SDset(iT,3),1))];
        TableOut{1+iT,4}  = num2str(round(p(iT),3));
        % TableOut{1+iT,4}  = [num2str(round(MeanSet(iT,1),1)) ' ? ' num2str(round(SDset(iT,1),1))];
    end
else
    for iT=1:length(CurrName)
        TableOut{iT,1}  = CurrName{iT};
        TableOut{iT,2}  = [num2str(round(MeanSet(iT,1),1)) ' ? ' num2str(round(SDset(iT,1),1))];
    end
end

