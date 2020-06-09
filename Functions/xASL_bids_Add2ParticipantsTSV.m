function xASL_bids_Add2ParticipantsTSV(DataIn, DataName, x, bOverwrite)
%xASL_bids_Add2ParticipantsTSV Add data to participants.tsv
%
% FORMAT: xASL_bids_Add2ParticipantsTSV(DataIn, DataName, x, bOverwrite)
%
% INPUT:
%   DataIn      - cell array with data to be added to participants.tsv, (REQUIRED)
%                 preferably with subjects/sessions (runs) as first & aecond columns, but
%                 this is not required. Number of datapoints should equal nSubjectsSessions.
%   DataName    - Name of the added variable/parameter (REQUIRED)
%   x           - struct containing the ExploreASL environment parameters
%                 (REQUIRED)
%   bOverwrite  - boolean specifying to overwrite pre-existing columns with
%                 the same header (OPTIONAL, DEFAULT = true) 
%                 (pre-existing participants.tsv is always overwritten)
%
% OUTPUT: n/a
% OUTPUT FILE:  /MyStudy/participants.tsv
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function adds metadata/statistical variables to the
% participants.tsv in the root/analysis folder, by the following steps.
% Note that this function will fail if the nSubjectsSessions is not the
% same for ExploreASL/data export, and participants.tsv!
%
% 1) Admin - Validate that there are not too many columns
% 2) Admin - Detect nSubjectsSessions
% 3) Admin - Load pre-existing participants.tsv or create one
% 4) Validate that subjects/session columns are equal with those in CellArray
% 5) Remove any existing columns
% 6) Add data to CellArray
% 7) Write data to participants.tsv
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_bids_Add2ParticipantsTSV(MeanMotion, 'MeanMotion', x);
% __________________________________
% Copyright 2015-2020 ExploreASL


PathTSV = fullfile(x.D.ROOT, 'participants.tsv');


%% 1) Admin - Validate that there are not too many columns
if nargin<4 || isempty(bOverwrite)
    bOverwrite = true;
end

if size(DataIn, 2)>3
    warning('Too many or few columns, skipping');
end

%% 2) Admin - Detect nSubjectsSessions
if size(DataIn, 1)==x.nSubjectsSessions && size(DataIn, 1)~=x.nSubjects
    HasSessions = true;
elseif size(DataIn, 1)==x.nSubjects
    HasSessions = false;
else
    warning('Incorrect nDatapoints, needs to be either nSubjects or nSubjectsSessions, skipping');
    return;
end

%% 3) Admin - Load pre-existing participants.tsv or create one
if exist(PathTSV, 'file')
    % 3A) Load & sort columns from participants.tsv
    CellArrayOrig = xASL_tsvRead(PathTSV);
    SubjectIndex = find(cellfun(@(y) ~isempty(regexp(y,'^(participant|subject).*id*$')), lower(CellArrayOrig(1,:))));
    SessionIndex = find(cellfun(@(y) ~isempty(regexp(y,'^(session).*(id)*$')), lower(CellArrayOrig(1,:))));
    
    % 3B) Check if participants.tsv has correct number of subjects/sessions
    if size(CellArrayOrig, 1)-1~=x.nSubjectsSessions
        warning([PathTSV ' doesnt have correct nSubjectsSessions, skipping']);
        return;
    end
    
    if isempty(SubjectIndex)
        warning('Missing participant_id column in pre-existing participant.tsv');
        return;
    elseif isempty(SessionIndex)
        warning('Missing session_id column in pre-existing participant.tsv');
        return;
    else % Sort columns to start with participant_id & session_id
        CellArray = CellArrayOrig(:, [SubjectIndex SessionIndex]);
        % Remove subjects & sessions from original CellArray
        CellArrayOrig = CellArrayOrig(:, [1:SubjectIndex-1 SubjectIndex+1:end]);
        SessionIndex = find(cellfun(@(y) ~isempty(regexp(y,'^(session).*(id)*$')), lower(CellArrayOrig(1,:))));
        CellArrayOrig = CellArrayOrig(:, [1:SessionIndex-1 SessionIndex+1:end]);
        CellArray = [CellArray CellArrayOrig];
    end
    
else
    % 3C) Create required columns subjects & sessions from x.SUBJECTS & x.SESSIONS

    for iSession=1:x.nSessions
        CellArray(1+iSession:x.nSessions:x.nSubjectsSessions+1, 1) = x.SUBJECTS;
    end
    CellArray(2:x.nSubjectsSessions+1, 2) = repmat(x.SESSIONS,[x.nSubjects, 1]);
end

% Force these names
CellArray{1,1} = 'participant_id';
CellArray{1,2} = 'session';

SubjectIndex = 1;
SessionIndex = 2;

%% -------------------------------------------------------------------------------------------
%% 4) Validate that subjects/session columns are equal with those in CellArray
if size(DataIn, 2)>1  % If to be included data have no sessions

    % Specify subjectIDs for metadata to be merged
    SubjectsShouldBe = DataIn(:,1);
    if size(DataIn, 1)==x.nSubjectsSessions
        SubjectsAre = CellArray(2:end, SubjectIndex);
    elseif size(DataIn, 1)==x.nSubjects
        SubjectsAre = CellArray(2:x.nSessions:end, SubjectIndex);
    end
        
    % Compare these for each subject
    for iSubject=1:length(SubjectsAre)
        if ~strcmp(SubjectsAre{iSubject}, SubjectsShouldBe{iSubject})
            warning('Invalid subject column in DataIn, skipping');
            return;
        end
    end
    
elseif size(DataIn, 2)>2 % If to be included data have sessions
    % Specify sessionIDs for metadata to be merged
    SessionsShouldBe = DataIn(:,2);
    SessionsAre = CellArray(2:end, SessionIndex);
    
    % Compare these for each session
    for iSession=1:length(SessionsAre)
        if ~strcmp(SessionsAre{iSession}, SessionsShouldBe{iSession})
            warning('Invalid session column in DataIn, skipping');
            return;
        end
    end
end

%% -------------------------------------------------------------------------------------------
%% 5) Remove any existing columns
% Find columns with the same header
HeaderIndices = strcmp(CellArray(1,:), DataName);
bAlreadyExisted = sum(HeaderIndices)>0;
if bAlreadyExisted && ~bOverwrite
    warning('Data column already existed, skipping');
    return;
elseif bAlreadyExisted && bOverwrite
    % Remove any pre-existing data with the same header
    CellArray = CellArray(:, ~HeaderIndices);
end

%% -------------------------------------------------------------------------------------------
%% 6) Add data to CellArray
CellArray{1, end+1} = DataName;
if size(DataIn, 1)==x.nSubjectsSessions % add session/run-specific data
    CellArray(2:end,end) = DataIn(:,end);
elseif size(DataIn, 1)==x.nSubjects % copy subject-specific data over multiple sessions/runs
    for iSession=1:x.nSessions
        CellArray(1+iSession:x.nSessions:x.nSubjectsSessions+1, end) = DataIn(:,end);
    end
else
    warning('Incorrect nSubjects/nSessions size of DataIn');
end

%% -------------------------------------------------------------------------------------------
%% 7) Write data to participants.tsv
xASL_tsvWrite(CellArray, PathTSV, 1);


end