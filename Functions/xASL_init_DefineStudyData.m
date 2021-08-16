function [x] = xASL_init_DefineStudyData(x)
%xASL_init_DefineStudyData Define study subjects/parameters for this
% pipeline run.
%
% FORMAT: [x] = xASL_init_DefineStudyData(x)
%
% INPUT:
%   x   - struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x   - struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This initialization wrapper initializes the parameters for
%               this pipeline run, i.e. subjects, sessions (runs), timepoints (visits),
%               exclusions, sites, cohorts etc.
%
%               Note that ASL sessions are defined here as what BIDS calls "runs".
%
%               The "longitudinal_Registration functions here manage different
%               TimePoints, which is what BIDS calls "visits".
%               With different structural scans, from the same participant. This is
%               managed by subject name suffixes _1 _2 _n, and can be used for comparing
%               visits in the population module, or running SPM's longitudinal within-subject
%               registration.
%
%               Parallelization is allowed here by calling ExploreASL different times,
%               where it divides the subjects/images for processing across the nWorkers,
%               using iWorker as the reference for the part that the current ExploreASL
%               call will process. This requires having a Matlab license that can be
%               started multiple times on a server, or alternatively running the
%               ExploreASL compilation, and doesn't require the Matlab parallel toolbox.
%
%               This function exists from the following parts:
%
% 1. Manage included & excluded subjects
% 2. Create dummy defaults (exclusion list, ASL sessions)
% 3. Create list of total baseline & follow-up subjects, before exclusions
% 4. Create TimePoint data-lists
% 5. Manage exclusions
% 6. Add sessions as statistical variable, if they exist
% 7. Parallelization: If running parallel, select cases for this worker
% 8. Add Longitudinal TimePoints (different T1 volumes) as covariate/set, after excluding
% 9. Load & add statistical variables
% 10. Make x.S.SetsOptions horizontal if vertical by transposing
% 11. Create single site dummy, if there were no sites specified
% 12. Check for time points sets
% 13. Create list of baseline & follow-up subjects (i.e. after exclusion)
% 14. Check what excluded from which TimePoints
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE: x = xASL_init_DefineStudyData(x);
% __________________________________
% Copyright 2015-2020 ExploreASL


% ------------------------------------------------------------------------------------------------
%% Admin
if ~isfield(x,'dataset')
    x.dataset = struct; 
end
if ~isfield(x.dataset,'name')
    x.dataset.name = ''; 
end


% ------------------------------------------------------------------------------------------------
%% 1) Manage included & excluded subjects
% Get a list of subject ID's by querying the folder names, excluding exclusions

if ~isfield(x.dataset,'subjectRegexp')
    warning('No subject regular expression detected');
    fprintf('Check if the correct DataPar file was loaded\n');
end

if isfield(x.dataset,'ForceInclusionList')
    % This is an option if you want to select subjects yourself,
    % instead of using all the subjects that comply with the regular expression
    x.dataset.TotalSubjects = x.dataset.ForceInclusionList';
else
    if ~isfield(x.dataset,'subjectRegexp')
        error('Please define the subjectRegexp of your dataset...');
    end
    % First escape double escaping
    x.dataset.subjectRegexp = strrep(x.dataset.subjectRegexp,'\\','\');
    % add the visit-postfix as option
    x.dataset.subjectRegexp = strrep(x.dataset.subjectRegexp,'$','');
    x.dataset.subjectRegexp = [x.dataset.subjectRegexp '(|_\d*)$'];
    % Then load subjects
    x.dataset.TotalSubjects = sort(xASL_adm_GetFileList(x.D.ROOT, x.dataset.subjectRegexp, 'List', [0 Inf], true)); % find dirs
end

x.dataset.nTotalSubjects = length(x.dataset.TotalSubjects);

% ------------------------------------------------------------------------------------------------
%% 2) Create dummy defaults (exclusion list, ASL sessions)
fprintf('Automatically defining sessions...\n');
if ~isfield(x,'dataset')
    x.dataset = struct;
end
if ~isfield(x.dataset,'exclusion') % default no exclusions
    x.dataset.exclusion = {''};
end

if isfield(x,'SESSIONS') && isstruct(x.SESSIONS)
    warning('Invalid x.SESSIONS structure, replacing this now');
    fprintf('%s\n', 'Check that the correct number of sessions were processed');
    x = rmfield(x,'SESSIONS');
end

if ~isfield(x,'SESSIONS')
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

% ------------------------------------------------------------------------------------------------
%% 3) Create list of total baseline & follow-up subjects, before exclusions
x.dataset.nSessions = length(x.SESSIONS);
x.SUBJECTS = x.dataset.TotalSubjects; % temporarily for xASL_init_LongitudinalRegistration
x.nSubjects = length(x.SUBJECTS);

if isempty(x.SUBJECTS)
    fprintf('No subjects found\n');
    fprintf('Please check the sub_regexp in your Data Parameter File\n');
    fprintf('This should match with the subject folders inside the ROOT folder\n');
    fprintf(['This was ' x.D.ROOT '\n']);
    error('No subjects defined, x.SUBJECTS was empty');
end

[~, TimePoint] = xASL_init_LongitudinalRegistration(x);

% ------------------------------------------------------------------------------------------------
%% 4) Create TimePoint data-lists
for iT=unique(TimePoint)'
    x.dataset.TimePointTotalSubjects{iT} = '';
end

for iSubj=1:x.dataset.nTotalSubjects
    iSess=1;
    iSubjSess = (iSubj-1)*x.dataset.nSessions + iSess;
    x.dataset.TimePointTotalSubjects{TimePoint(iSubjSess)}{end+1} = x.dataset.TotalSubjects{iSubj};
end
    
x = rmfield(x,'SUBJECTS');
x = rmfield(x,'nSubjects');

% ------------------------------------------------------------------------------------------------
%% 5) Manage exclusions
if     ~isfield(x.dataset,'exclusion')
        x.dataset.exclusion = {''};
        nExclusion = 0;
elseif  isempty(x.dataset.exclusion)
        x.dataset.exclusion = {''};    
        nExclusion = 0;
else
        nExclusion = length(x.dataset.exclusion);
end

if ~iscell(x.dataset.exclusion)
    x.dataset.exclusion = {x.dataset.exclusion};
end

if  nExclusion==1 && isempty(x.dataset.exclusion{1})
    nExclusion = 0;
end

x.dataset.ExcludedSubjects = '';
x.SUBJECTS = '';
for iSubject=1:x.dataset.nTotalSubjects
    excl=0;
    for j=1:nExclusion % Check if subject should be excluded
        if  strcmp(x.dataset.TotalSubjects{iSubject},x.dataset.exclusion{j})
            excl=1;
            x.dataset.ExcludedSubjects{end+1} = x.dataset.TotalSubjects{iSubject};
        end
    end
    
    if ~isempty(regexp(x.dataset.TotalSubjects{iSubject},'^(dartel|lock|Population)$'))
         % This is not a subject but a pipeline folder
         ListNoPipelineDir(iSubject) = 0;
    elseif ~excl
         % include subject if not to be excluded, and if it doesn't have
         % same name as pipeline directories     
         x.SUBJECTS{end+1}  = x.dataset.TotalSubjects{iSubject};
         x.dataset.TotalInclusionList( (iSubject-1)*x.dataset.nSessions+1:iSubject*x.dataset.nSessions,1) = 1;
         ListNoPipelineDir(iSubject) = 1;
    else
        x.dataset.TotalInclusionList( (iSubject-1)*x.dataset.nSessions+1:iSubject*x.dataset.nSessions,1) = 0;
        ListNoPipelineDir(iSubject) = 1;
    end
end

% Remove pipeline dirs from list
if sum(ListNoPipelineDir)>0 % if there are any subjects
    x.dataset.TotalSubjects = x.dataset.TotalSubjects(logical(ListNoPipelineDir));
end

if ~isfield(x,'SUBJECTS')
    fprintf('%s\n','No subjects found');
    x.SUBJECTS = [];
end

x.nSubjects = length(x.SUBJECTS);
x.dataset.nTotalSubjects = length(x.dataset.TotalSubjects);
x.dataset.nExcluded = x.dataset.nTotalSubjects - x.nSubjects;
x.dataset.nSessions = length( x.SESSIONS );
x.dataset.nSubjectsSessions = x.nSubjects .* x.dataset.nSessions;

x.dataset.nTimePointsTotal = length(x.dataset.TimePointTotalSubjects);
for iT=1:x.dataset.nTimePointsTotal
    x.dataset.nTimePointTotalSubjects(iT) = length(x.dataset.TimePointTotalSubjects{iT});
end

% ------------------------------------------------------------------------------------------------
%% 6) Add sessions as statistical variable

% Predefine SETS to avoid empty SETS & import predefined session settings as set settings
x.S.SetsName{1} = 'session';
x.S.SetsOptions = '';
x.S.SetsID = 0;
x.S.Sets1_2Sample(1) = 1; % sessions are always paired observations (e.g. belonging to same individual, looking at intra-individual changes)

% NOT SURE IF WE STILL USE THIS, WE MAY PHASE THIS OUT
% Define SetsOptions for sessions
if ~isfield(x,'session')
    x.session = '';
end

if  isfield(x.session,'options')
    x.S.SetsOptions{1} = x.session.options; % with session options (e.g. morning, evening 2nd morning)
else
    for iS=1:x.dataset.nSessions
        x.S.SetsOptions{1}{iS} = ['ASL_' num2str(iS)];
    end
end        

% Create ID/data for session numbers
for iSubj=1:x.nSubjects
    for iSess=1:x.dataset.nSessions
        iSubjSess = (iSubj-1)*x.dataset.nSessions+iSess;

        x.S.SetsID(iSubjSess,1) = iSess;
    end
end

% ------------------------------------------------------------------------------------------------
%% 7) Parallelization: If running parallel, select cases for this worker
if x.opts.nWorkers>1
    nSubjPerWorker = ceil(x.nSubjects/x.opts.nWorkers); % ceil to make sure all subjects are processed
    nSubjSessPerWorker = nSubjPerWorker*x.dataset.nSessions;
    iStartSubject = (x.opts.iWorker-1)*nSubjPerWorker+1;
    iEndSubject = min(x.opts.iWorker*nSubjPerWorker, x.nSubjects);
    iStartSubjectSession = (x.opts.iWorker-1)*nSubjSessPerWorker+1;
    iEndSubjectSession = min(x.opts.iWorker*nSubjSessPerWorker, x.dataset.nSubjectsSessions);
    
    if iStartSubject>x.nSubjects
        fprintf('Closing down this worker, had too many workers');
        exit;
    end
    
    % Adapt SUBJECTS
    x.SUBJECTS = x.SUBJECTS(iStartSubject:iEndSubject);
    x.nSubjects = length(x.SUBJECTS);
    x.dataset.nSubjectsSessions = x.nSubjects*x.dataset.nSessions;

    % Adapt SETSID (covariants)
    if isfield(x.S,'SetsID') && ~isempty(x.S.SetsID)
        x.S.SetsID = x.S.SetsID(iStartSubjectSession:iEndSubjectSession, :);
    end
    
    fprintf(['I am worker ' num2str(x.opts.iWorker) '/' num2str(x.opts.nWorkers) '\n']);
    fprintf(['I will process subjects ' num2str(iStartSubject) '-' num2str(iEndSubject) '\n']);
end

% ------------------------------------------------------------------------------------------------
%% 8) Add Longitudinal TimePoints (different T1 volumes) as covariate/set, after excluding
% TimePoints should be indicated as different subjects, with a _1 _2 _3 _n suffix

[SubjectNList, TimePoint] = xASL_init_LongitudinalRegistration(x);

if isfield(x.S,'SetsID')
    iSetLong_TP = size(x.S.SetsID,2)+1;
else
    iSetLong_TP = 1;
end

x.S.SetsName{iSetLong_TP} = 'LongitudinalTimePoint';
x.S.SetsName{iSetLong_TP+1} = 'SubjectNList';
x.S.SetsID(:,iSetLong_TP) = TimePoint;
x.S.SetsID(:,iSetLong_TP+1) = SubjectNList;

for iT=1:max(TimePoint)
    x.S.SetsOptions{iSetLong_TP}{iT} = ['TimePoint_' num2str(iT)];
end

x.S.SetsOptions{iSetLong_TP+1} = {'SubjectNList'};
x.S.Sets1_2Sample(iSetLong_TP:iSetLong_TP+1) = 1; % Longitudinal TimePoints are always paired observations (e.g. belonging to same individual, looking at intra-individual changes)    
x.S.iSetLong_TP = iSetLong_TP;




% ================================================================================================
%% 9) Load & add statistical variables
% Keep some space here, to easily find this for debugging
x = xASL_init_LoadMetadata(x);











% ================================================================================================
%% 10) Make x.S.SetsOptions horizontal if vertical by transposing
if isfield(x.S,'SetsOptions')
    for iCell=1:length(x.S.SetsOptions)
        if size(x.S.SetsOptions{iCell},1)>size(x.S.SetsOptions{iCell},2)
            x.S.SetsOptions{iCell} = x.S.SetsOptions{iCell}';
        end
    end
else
    x.S.SetsOptions = {};
end


% ------------------------------------------------------------------------------------------------
%% 11) Create single site dummy, if there were no sites specified
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


% ------------------------------------------------------------------------------------------------
%% 12) Check for time points sets
% Search for a variable specifying longitudinal time points
if isfield(x.S, 'SetsName')
    % Search for a variable specifying age at scanning
    iS = find(strcmp(x.S.SetsName,'Age'));
    if ~isempty(iS)
        x.S.iSetAge=iS;
    end
end

% ------------------------------------------------------------------------------------------------
%% 13) Create list of baseline & follow-up subjects (i.e. after exclusion)
for iCell=1:length(x.dataset.TimePointTotalSubjects)
    x.dataset.TimePointSubjects{iCell} = '';
end
    
for iS=1:x.nSubjects
    iSession = 1; % append to accommodate sessions in SetsID
    iSubjSess = ((iS-1)*x.dataset.nSessions)+iSession;
    CurrentTimePoint = x.S.SetsID(iSubjSess,x.S.iSetLong_TP);
    if length(x.dataset.TimePointSubjects)<CurrentTimePoint % if this cell didn't exist yet
        x.dataset.TimePointSubjects{CurrentTimePoint} = ''; 
    end    
    x.dataset.TimePointSubjects{CurrentTimePoint}{end+1} = x.SUBJECTS{iS};
end
    
x.dataset.nTimePoints = length(x.dataset.TimePointSubjects);
for iT=1:x.dataset.nTimePoints
    if  length(x.dataset.TimePointSubjects)<iT
        % if an excluded later volume led to different composition
        % of TotalSubjects (i.e. before exclusion) & Subjects (i.e. after
        % exclusion)
        x.dataset.nTimePointSubjects(iT) = 0; 
    else
        x.dataset.nTimePointSubjects(iT) = length(x.dataset.TimePointSubjects{iT});
    end
end

% ------------------------------------------------------------------------------------------------
%% 14) Check what excluded from which TimePoints
for iT=1:x.dataset.nTimePoints
    x.dataset.TimePointExcluded{iT} = '';
end

for iE=1:x.dataset.nExcluded
    FoundE = 0;
    iT = 0;
    while ~FoundE % excluded subject not found in previous TimePoint
        iT=iT+1;  % go to next TimePoint
        
        iS  = find(strcmp(x.dataset.TimePointTotalSubjects{iT}, x.dataset.ExcludedSubjects{iE}));
        if ~isempty(iS)
            x.dataset.TimePointExcluded{iT}{end+1} = x.dataset.TimePointTotalSubjects{iT}{iS};
            FoundE = 1;
        end
    end
end

for iT=1:x.dataset.nTimePoints
    x.dataset.nTimePointExcluded(iT) = length(x.dataset.TimePointExcluded{iT});
end


end
