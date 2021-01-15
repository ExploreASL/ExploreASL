function [SubjectNlist, TimePoint, IsSubject, SubjectID_FirstVolume] = xASL_init_LongitudinalRegistration(x)
%xASL_init_LongitudinalRegistration Initialization of longitudinal registration
%
% FORMAT: [SubjectNlist, TimePoint, IsSubject, SubjectID_FirstVolume] = xASL_init_LongitudinalRegistration(x)
%
% INPUT:
%   x       - struct containing pipeline environment parameters (REQUIRED)
%
% OUTPUT:
%   SubjectNlist            - struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%   TimePoint               - TimePoint counter (1, 2, etc)
%   IsSubject               - timepoint index of current subject
%   SubjectID_FirstVolume   - name/ID of the (first visit/timepoint) of the subject
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function initializes the longitudinal registration for ExploreASL,
% which implements the SPM longitudinal registration.
%
% This function recognizes and defines visits (i.e. multiple scans per
% subject) in the initialization of ExploreASL, based on the suffixes _1 _2
% _n in the subject names (identified as foldernames).
%
% Specifically, this function is called in the registration modules LongReg and DARTEL,
% the first carrying out within-subject registration and 
% the second between-subject registration, based on the first time point
% only.
% For the first function, we specify here a list of visits/timepoints that
% should be registered longitudinally, for the second function we specify a
% list of first visits only, as the between-subject registration in
% ExploreASL is based on the first scan (as opposed to the average
% subject's scan).
%
% This function runs the following steps:
%
% 1. Get TimePoint-list (list of visits)
% 2. Find subject IDs
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [SubjectNlist, TimePoint, IsSubject, SubjectID_FirstVolume] = xASL_init_LongitudinalRegistration(x);
% __________________________________
% Copyright 2015-2020 ExploreASL


%% Administration
% if no sessions are defined yet, assume one session only
if ~isfield(x,'nSessions')
    x.nSessions = length(x.SESSIONS);
end

%%  1) Get TimePoint-list (list of visits)
% Including time point list with time points that exist in the data (e.g. if only volumes 3 & 4 exist for a subject, these will be called 1 & 2
% if it didn't exist yet, create it
% NB: volumes _3 _4 _5 will be called TimePoint 3 4 5 here

if isfield(x,'iSetLong_TP')
    TimePoint = x.S.SetsID(:,x.iSetLong_TP);
    SubjectNameList = x.SubjectNameList;
    SubjectNlist = x.S.SetsID(:,x.iSetLong_TP+1);
else
    for iS=1:x.nSubjects
        % First check whether if there is a TimePoint indication
        % With a maximum of 99 TimePoints, the underscore should be 
        % a maximum of 2 indices below the end of the SubjectName
        clear CurrentTimePoint
        Indices = find(x.SUBJECTS{iS}=='_');

        if ~isempty(Indices)
            if Indices(end)>=length(x.SUBJECTS{iS})-2
                % This is a TimePoint indication
                CurrentTimePoint = str2num( x.SUBJECTS{iS}(Indices(end)+1:end));
                CurrentSubjName  = x.SUBJECTS{iS}(1:Indices(end)-1);
            end
        end

        if ~exist('CurrentTimePoint','var') || isempty(CurrentTimePoint) || ~isnumeric(CurrentTimePoint)
            CurrentTimePoint = 1;
            CurrentSubjName = x.SUBJECTS{iS};
            % if no index suffix, then assume first TimePoint (default)
        end

        iSubjSessAll = [(iS-1)*x.nSessions+1:iS*x.nSessions];
        PreviousSubjSess =  (iS-1)*x.nSessions;
        SubjectNameList(iSubjSessAll,1) = {CurrentSubjName};

        if ~exist('SubjectNlist','var')
            % for first case
            SubjectNlist(iSubjSessAll,1) = 1;
            TimePoint(iSubjSessAll,1) = 1;
        elseif ~strcmp(CurrentSubjName,SubjectNameList{PreviousSubjSess,1})
            % if new subject, add to list as first time point
            TimePoint(iSubjSessAll,1) = 1; % first time point
            SubjectNlist(iSubjSessAll,1) = SubjectNlist(PreviousSubjSess,1) + 1; % previous subject + 1 == new value
        else
            % if same subject, add to list as new time point
            TimePoint(iSubjSessAll,1) = TimePoint(PreviousSubjSess,1) + 1; % previous subject + 1 == new time point
            SubjectNlist(iSubjSessAll,1) = SubjectNlist(PreviousSubjSess,1); % same as last time point
        end

        if ~isnumeric(CurrentTimePoint) || isempty(CurrentTimePoint) % sanity check
            error(['Subject-name ' x.SUBJECTS{iS} ' contains a underscore which is reserved for longitudinal TimePoints 1 2 3 etc']);
        end
    end
end

%% 2) Find subject IDs
% The pipeline should only perform longitudinal registration if this is the first time point
% Therefore, first find all volumes with the same subject name

if isfield(x, 'P')
    if isfield(x.P, 'SubjectID')
        iSubject = find(strcmp(x.SUBJECTS, x.P.SubjectID));
        % Subject session id number for the first session,, to accommodate that x.S.SetsID can contain sessions
        iSubjSess = (iSubject-1) * x.nSessions+1;  

        Indices = find(x.P.SubjectID=='_');
        if ~isempty(Indices)
            if Indices(end)>=length(x.P.SubjectID)-2
                Indices = Indices(end);
                SubjName = x.P.SubjectID(1:Indices-1);
            end
        else
        end

        if ~exist('SubjName','var')
            SubjName = x.P.SubjectID;
        end

        Indices = find(strcmp(SubjectNameList,SubjName));
        VolumeList = zeros(size(SubjectNameList,1),2);
        VolumeList(Indices,1:2) = [repmat(1,[length(Indices),1]) TimePoint(Indices)];

        % Then only perform longitudinal registration if our volume is the
        % lowest index, and if there are multiple indices    
        VolumeN = unique(sort(nonzeros(VolumeList(:,2))'));
        IsSubject = find(VolumeN==TimePoint(iSubjSess,1)); % give the volume index (e.g. 2, 3 etc)

        % find first volume/time point
        IndexFirstSubject = find(VolumeList(:,2)==min(VolumeN));
        iSubj = xASL_adm_ConvertSubjSess2Subj_Sess(x.nSessions, IndexFirstSubject);
        SubjectID_FirstVolume = x.SUBJECTS{iSubj};
    end
end
    
    
end