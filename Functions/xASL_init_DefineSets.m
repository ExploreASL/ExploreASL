function x = xASL_init_DefineSets(x)
%xASL_init_DefineSets Part of master script ExploreASL, defining sets
% from sessions and pre-defined groups
% TotalSubjects == subjects before processing defined exclusions
% ExploreASL, HJMM Mutsaerts, 2016
    
%% Manage included & excluded subjects
% Get a list of subject ID's by querying the folder names, excluding exclusions

if  isfield(x,'ForceInclusionList')
    % This is an option if you want to select subjects
    x.TotalSubjects = x.ForceInclusionList';
else
    x.TotalSubjects = sort(xASL_adm_GetFileList(x.D.ROOT, x.subject_regexp, 'List', [0 Inf], true)); % find dirs
end


x.nTotalSubjects = length(x.TotalSubjects);

%% Create dummy defaults
if ~isfield(x,'exclusion') % default no exclusions
    x.exclusion = {''};
end

if ~isfield(x,'SESSIONS') % can also be defined in DataPar.mat
    x.SESSIONS = '';
    SessionPathList = xASL_adm_GetFileList(x.D.ROOT, '^(ASL|func)_\d*$', 'FPListRec', [0 Inf],1);
    for iSess=1:length(SessionPathList)
        [~, x.SESSIONS{end+1}]  = fileparts(SessionPathList{iSess});
    end
    x.SESSIONS = unique(x.SESSIONS);
    if  isempty(x.SESSIONS)
        x.SESSIONS{1} = 'ASL_1'; % default session
    end
end

%% Create list of total baseline & follow-up subjects, before exclusions
x.nSessions = length(x.SESSIONS);
x.SUBJECTS = x.TotalSubjects; % temporarily for xASL_init_LongitudinalRegistration
x.nSubjects = length(x.SUBJECTS);

if  isempty(x.SUBJECTS)
    error('No subjects defined, x.SUBJECTS was empty');
end

[~, ~, TimePoint] = xASL_init_LongitudinalRegistration( x );


%% Create TimePoint data-lists
for iT=unique(TimePoint)'
    x.TimePointTotalSubjects{iT} = '';
end

for iSubj=1:x.nTotalSubjects
    iSess=1;
    iSubjSess = (iSubj-1)*x.nSessions + iSess;
    x.TimePointTotalSubjects{TimePoint(iSubjSess)}{end+1} = x.TotalSubjects{iSubj};
end
    
x = rmfield(x,'SUBJECTS');
x = rmfield(x,'nSubjects');

%% Manage exclusions

if     ~isfield(x,'exclusion')
        x.exclusion = {''};
        nExclusion = 0;
elseif  isempty(x.exclusion)
        x.exclusion = {''};    
        nExclusion = 0;
else
        nExclusion = length(x.exclusion);
end

if ~iscell(x.exclusion)
    x.exclusion = {x.exclusion};
end

if  nExclusion==1 && isempty(x.exclusion{1})
    nExclusion = 0;
end

x.ExcludedSubjects = '';
x.SUBJECTS = '';
for iSubject=1:x.nTotalSubjects
    excl=0;
    for j=1:nExclusion % Check if subject should be excluded
        if  strcmp(x.TotalSubjects{iSubject},x.exclusion{j})
            excl=1;
            x.ExcludedSubjects{end+1} = x.TotalSubjects{iSubject};
        end
    end
    
    if      strcmp(x.TotalSubjects{iSubject},'dartel') || strcmp(x.TotalSubjects{iSubject},'lock')
            % Remove this pseudo-subject from the list
            ListNoPipelineDir(iSubject) = 0;
            
    elseif ~excl
            % include subject if not to be excluded, and if it doesn't have
            % same name as pipeline directories
            x.SUBJECTS{end+1}  = x.TotalSubjects{iSubject};
            x.TotalInclusionList( (iSubject-1)*x.nSessions+1:iSubject*x.nSessions,1)=1;
            ListNoPipelineDir(iSubject) = 1;
    else
            x.TotalInclusionList( (iSubject-1)*x.nSessions+1:iSubject*x.nSessions,1)=0;
            ListNoPipelineDir(iSubject) = 1;
    end
end

% Remove pipeline dirs from list
NextN = 1;
if  sum(ListNoPipelineDir)>0
    TotalSubjectsTemp = x.TotalSubjects;
    x = rmfield(x,'TotalSubjects');
    for iS=1:length(TotalSubjectsTemp)
        if  ListNoPipelineDir(iS)
            x.TotalSubjects{NextN} = TotalSubjectsTemp{iS};
            NextN = NextN+1;
        end
    end
end
           

x.SUBJECTS = sort(x.SUBJECTS);

clear iSubject j

if ~isfield(x,'SUBJECTS')
    fprintf('%s\n','No subjects found');
    x.SUBJECTS = [];
end

x.nSubjects = length(x.SUBJECTS);
x.nTotalSubjects = length(x.TotalSubjects);
x.nExcluded = x.nTotalSubjects - x.nSubjects;
x.nSessions = length( x.SESSIONS );
x.nSubjectsSessions = x.nSubjects .* x.nSessions;

x.nTimePointsTotal = length(x.TimePointTotalSubjects);
for iT=1:x.nTimePointsTotal
    x.nTimePointTotalSubjects(iT) = length(x.TimePointTotalSubjects{iT});
end

%% Add sessions as statistical variable, if they exist
if  x.nSessions>1  % if there are sessions (more than 1 session), then sessions=1st set

    % Predefine SETS to avoid empty SETS & import predefined session settings as set settings
    x.S.SetsName{1} = 'session';
    x.S.SetsOptions = '';
    x.S.SetsID = 0;
    x.S.Sets1_2Sample(1) = 1; % sessions are always paired observations (e.g. belonging to same individual, looking at intra-individual changes)
    
    % Define SetsOptions for sessions
    if ~isfield(x,'session')
        x.session = '';
    end
        
    if  isfield(x.session,'options')
        x.S.SetsOptions{1} = x.session.options; % with session options (e.g. morning, evening 2nd morning)
    else
        for iS=1:x.nSessions
            x.S.SetsOptions{1}{iS} = ['Session_' num2str(iS)];
        end
    end        
    

    
    %% Create ID/data for session numbers
    for iSubj=1:x.nSubjects
        for iSess=1:x.nSessions
            iSubjSess = (iSubj-1)*x.nSessions+iSess;
            
            x.S.SetsID(iSubjSess,1) = iSess;
        end
    end
end





%% PARALLELIZATION: If running parallel, select cases for this worker
if x.nWorkers>1
    nSubjPerWorker = ceil(x.nSubjects/x.nWorkers); % ceil to make sure all subjects are processed
    iStart = (x.iWorker-1)*nSubjPerWorker+1;
    iEnd = min(x.iWorker*nSubjPerWorker, x.nSubjects);
    
    if iStart>x.nSubjects
        fprintf('Closing down this worker, had too many workers');
        exit;
    end
    
    % Adapt SUBJECTS
    x.SUBJECTS = x.SUBJECTS(iStart:iEnd);
    x.nSubjects = length(x.SUBJECTS);
    x.nSubjectsSessions = x.nSubjects*x.nSessions;
    
    % Adapt SETSID (covariants)
    if isfield(x.S,'SetsID') && ~isempty(x.S.SetsID)
        x.S.SetsID = x.S.SetsID(iStart:iEnd, :);
    end
    
    fprintf(['I am worker ' num2str(x.iWorker) '/' num2str(x.nWorkers) '\n']);
    fprintf(['I will process subjects ' num2str(iStart) '-' num2str(iEnd) '\n']);
end





%% 1) Add Longitudinal TimePoints (different T1 volumes) as covariate/set, AFTER EXCLUSION
% TimePoints should be indicated as different subjects, with a _1 _2 _3 _n
% suffix

[SubjectNList, ~, TimePoint] = xASL_init_LongitudinalRegistration( x );



if isfield(x.S,'SetsID')
    iSetLong_TP = size(x.S.SetsID,2)+1;
else
    iSetLong_TP = 1;
end

x.S.SetsName{iSetLong_TP} = 'LongitudinalTimePoint';
x.S.SetsName{iSetLong_TP+1} = 'SubjectNList';
x.S.SetsID(:,iSetLong_TP) = TimePoint;
x.S.SetsID(:,iSetLong_TP+1) = SubjectNList;

for iT=1:max(TimePoint); x.S.SetsOptions{iSetLong_TP}{iT} = ['TimePoint_' num2str(iT)]; end

x.S.SetsOptions{iSetLong_TP+1} = {'SubjectNList'};
x.S.Sets1_2Sample(iSetLong_TP:iSetLong_TP+1) = 1; % Longitudinal TimePoints are always paired observations (e.g. belonging to same individual, looking at intra-individual changes)    
x.S.iSetLong_TP = iSetLong_TP;








%% LOAD STATS: Add statistical variables
x = xASL_init_LoadStatsData(x);








%% 2) Create single site dummy, if there were no sites specified
if isfield(x.S, 'SetsName')
    % Search for a Site variable
    for iS=1:length(x.S.SetsName)
        if  strcmp(x.S.SetsName{iS},'Site')
            x.S.iSetSite=iS;
        end
    end
end
if ~isfield(x.S,'iSetSite')
    x.S.SetsName{end+1} = 'Site';
    x.S.SetsOptions{end+1} = {'SingleSite'};
    x.S.SetsID(:,end+1) = 1;
    x.S.Sets1_2Sample(end+1) = 2;
    x.S.iSetSite = length(x.S.SetsName);
end   



%% Check for time points sets
% Search for a variable specifying longitudinal time points
if isfield(x.S, 'SetsName')
    % Search for a variable specifying age at scanning
    iS = find(strcmp(x.S.SetsName,'Age'));
    if ~isempty(iS)
        x.S.iSetAge=iS;
    end
end

%% Create list of baseline & follow-up subjects (i.e. after exclusion)
x.TimePointSubjects{1} = '';
    
for iS=1:x.nSubjects
    iSession = 1; % append to accommodate sessions in SetsID
    iSubjSess = ((iS-1)*x.nSessions)+iSession;
    CurrentTimePoint = x.S.SetsID(iSubjSess,x.S.iSetLong_TP);
    if length(x.TimePointSubjects)<CurrentTimePoint % if this cell didn't exist yet
        x.TimePointSubjects{CurrentTimePoint} = ''; 
    end    
    x.TimePointSubjects{CurrentTimePoint}{end+1} = x.SUBJECTS{iS};
end
    
x.nTimePoints = length(x.TimePointSubjects);
for iT=1:x.nTimePoints
    if  length(x.TimePointSubjects)<iT
        % if an excluded later volume led to different composition
        % of TotalSubjects (i.e. before exclusion) & Subjects (i.e. after
        % exclusion)
        x.nTimePointSubjects(iT) = 0; 
    else
        x.nTimePointSubjects(iT) = length(x.TimePointSubjects{iT});
    end
end


%% Check what excluded from which TimePoints
for iT=1:x.nTimePoints
    x.TimePointExcluded{iT} = '';
end

for iE=1:x.nExcluded
    FoundE = 0;
    iT = 0;
    while ~FoundE % excluded subject not found in previous TimePoint
        iT=iT+1;  % go to next TimePoint
        
        iS  = find(strcmp( x.TimePointTotalSubjects{iT}, x.ExcludedSubjects{iE} ));
        if ~isempty(iS)
            x.TimePointExcluded{iT}{end+1} = x.TimePointTotalSubjects{iT}{iS};
            FoundE = 1;
        end
    end
end

for iT=1:x.nTimePoints
    x.nTimePointExcluded(iT) = length(x.TimePointExcluded{iT});
end


end
