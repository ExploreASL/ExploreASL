function [SubjectNlist, SubjectNameList, TimePointCorr, IsSubject, VolumeList, VolumeN, SubjectID_FirstVolume ] = xASL_init_LongitudinalRegistration( x )
%xASL_init_LongitudinalRegistration Initialization of longitudinal registration
% To check whether or not we will run longitudinal registration
% IsVolume == 0 -> no Longitudinal Registration
% IsVolume == 1 -> do Longitudinal Registration, this is first TimePoint
% IsVolume  > 1 -> do Longitudinal Registration, this is a later TimePoint
% VolumeN provides the numbers of the TimePoints (this can be [1 2 3] for
% three subsequent TimePoints, but could also be [2 3] if the first
% TimePoint was excluded from ExploreASL

%% Administration
% if no sessions are defined yet, assume first session, e.g. when this
% function is used before sessions are defined in initialization
% ExploreASL
if ~isfield(x,'nSessions')
    x.nSessions   = length(x.SESSIONS);
end

%%  Get TimePoint-list, including time point list with time points that exist in the data (e.g. if only volumes 3 & 4 exist for a subject, these will be called 1 & 2
    % if it didn't exist yet, create it
    % NB: volumes _3 _4 _5 will be called TimePoint 3 4 5 here

if  isfield(x,'iSetLong_TP')
    TimePointCorr               = x.S.SetsID(:,x.iSetLong_TP);
    SubjectNameList             = x.SubjectNameList;
    SubjectNlist                = x.S.SetsID(:,x.iSetLong_TP+1);
else
    for iS=1:x.nSubjects
        % First check whether if there is a TimePoint indication
        % With a maximum of 99 TimePoints, the underscore should be 
        % a maximum of 2 indices below the end of the SubjectName
        clear CurrentTimePoint
        Indices     = find(x.SUBJECTS{iS}=='_');

        if ~isempty(Indices)
            if  Indices(end)>=length(x.SUBJECTS{iS})-2
                % This is a TimePoint indication
                CurrentTimePoint = str2num( x.SUBJECTS{iS}(Indices(end)+1:end) );
                CurrentSubjName  = x.SUBJECTS{iS}(1:Indices(end)-1);
            end
        end

        if ~exist('CurrentTimePoint','var')
            CurrentTimePoint    = 1;
            CurrentSubjName     = x.SUBJECTS{iS};
            % if no index suffix, then assume first TimePoint (default)
        end

        iSubjSessAll                    = [(iS-1)*x.nSessions+1:iS*x.nSessions];
        PreviousSubjSess                =  (iS-1)*x.nSessions;
        TimePoint(iSubjSessAll,1)       = CurrentTimePoint; % the actual value in the name
        SubjectNameList(iSubjSessAll,1) = {CurrentSubjName};

        if ~exist('SubjectNlist','var')
            % for first case
            SubjectNlist(iSubjSessAll,1)           = 1;
            TimePointCorr(iSubjSessAll,1)       = 1;
        elseif  ~strcmp(CurrentSubjName,SubjectNameList{PreviousSubjSess,1})
            % if new subject, add to list as first time point
            TimePointCorr(iSubjSessAll,1)  = 1; % first time point
            SubjectNlist(iSubjSessAll,1)      = SubjectNlist( PreviousSubjSess,1)+1; % previous subject + 1 == new value
        else
            % if same subject, add to list as new time point
            TimePointCorr(iSubjSessAll,1)  = TimePointCorr( PreviousSubjSess,1) + 1; % previous subject + 1 == new time point
            SubjectNlist(iSubjSessAll,1)      = SubjectNlist(PreviousSubjSess,1); % same as last time point
        end

        if ~isnumeric(CurrentTimePoint) || isempty(CurrentTimePoint) % sanity check
            error(['Subject-name ' x.SUBJECTS{iS} ' contains a underline which is reserved for longitudinal TimePoints 1 2 3 etc']);
        end
    end
end

%% Determine whether we will run longitudinal registration, get TimePoint numbers
% Find subject id number

% only perform longitudinal registration if this is the first time point
% Therefore, first find all volumes with the same subject name

if  isfield(x,'P')
    if  isfield(x.P,'SubjectID')
        iSubject    = find(strcmp(x.SUBJECTS,x.P.SubjectID));
        % Subject session id number for the first session,, to accommodate that x.S.SetsID can contain sessions
        iSubjSess   = (iSubject-1)*x.nSessions+1;  

        Indices     = find(x.P.SubjectID=='_');
        if ~isempty(Indices)
            if  Indices(end)>=length(x.P.SubjectID)-2
                Indices         = Indices(end);
                SubjName        = x.P.SubjectID(1:Indices-1);
            end
        else
        end

        if ~exist('SubjName','var')
            SubjName            = x.P.SubjectID;
        end


        Indices                     = find(strcmp(SubjectNameList,SubjName));
        VolumeList                  = zeros(size(SubjectNameList,1),2);
        VolumeList(Indices,1:2)     = [repmat(1,[length(Indices),1]) TimePointCorr(Indices)];

        % Then only perform longitudinal registration if our volume is the
        % lowest index, and if there are multiple indices    
        VolumeN                 = unique(sort(nonzeros(VolumeList(:,2))'));
        IsSubject               = find( VolumeN==TimePointCorr(iSubjSess,1)); % give the volume index (e.g. 2, 3 etc)

        % find first volume/time point
        IndexFirstSubject = find(VolumeList(:,2)==min(VolumeN));
        iSubj = xASL_adm_ConvertSubjSess2Subj_Sess(x.nSessions, IndexFirstSubject);
        SubjectID_FirstVolume   = x.SUBJECTS{iSubj};
    end
end
    
    
end

