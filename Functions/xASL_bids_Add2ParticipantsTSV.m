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
%   bOverwrite  - boolean specifying to overwrite pre-existing data values for subject/session with
%                 the same header (OPTIONAL, DEFAULT = true) 
%                 (pre-existing participants.tsv is always overwritten)
%
% OUTPUT: n/a
% OUTPUT FILE:  /MyStudy/participants.tsv
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function adds metadata/statistical variables to the
%               participants.tsv in the root/analysis folder, by the following steps.
%               This function will iterate over Data provided at DataIn and fill them
%               in the participants.tsv, overwriting if allowed.
%               Empty data is filled in as 'n/a', and the first column "participants_id"
%               is sorted for participants.
%
% This function runs the following steps:
%
% 1. Admin - Validate that there are not too many columns
% 2. Admin - Detect nSubjectsSessions
% 3. Admin - Load pre-existing participants.tsv or create one
% 4. Admin - Get column number of data
% 5. Add data to CellArray
% 6. Sort rows on subjects
% 7. Fill empty cells
% 8. Write data to participants.tsv
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
% This part checks if the DataIn has an invalid size.
% After this part, the columns should be 1) subject, 2) session/run, 3)
% data value
HadSessions = false;
if size(DataIn, 1)==x.dataset.nSubjectsSessions && size(DataIn, 1)~=x.nSubjects
    % A session column should exist
    if size(DataIn,2)<3
        warning('Session column missing, too few columns, skipping');
        return;
    else
        TempSessionNames = unique(DataIn(:,2));
        if ~max(strcmp(TempSessionNames, 'ASL_1'))
            warning('First session missing, might go wrong');
        end
        HadSessions = true;
    end
elseif size(DataIn, 1)==x.nSubjects
    % Add session column
    if size(DataIn, 2)==2
        DataIn(:,3) = DataIn(:,2);
        DataIn(:,2) = repmat({'ASL_1'}, [x.nSubjects 1]);
    elseif size(DataIn, 2)~=3
        warning('Too many or few columns, something went wrong, skipping');
        return;
    end
else
    warning('Incorrect nDatapoints, needs to be either nSubjects or nSubjectsSessions, skipping');
    return;
end

%% 3) Admin - Load pre-existing participants.tsv or create one
% This part ensures that a subject & session column exists.
% If a session column doesn't exist, it will be created.
if exist(PathTSV, 'file')
    % 3A) Load & sort columns from participants.tsv
    CellArrayOrig = xASL_tsvRead(PathTSV);
    SubjectIndex = find(cellfun(@(y) ~isempty(regexp(y,'^(participant|subject).*id*$')), lower(CellArrayOrig(1,:))));
    SessionIndex = find(cellfun(@(y) ~isempty(regexp(y,'^(session).*(id)*$')), lower(CellArrayOrig(1,:))));
    
    % 3B) Check participants.tsv columns
    if isempty(SubjectIndex)
        warning('Missing participant_id column in pre-existing participant.tsv');
        return;
    elseif isempty(SessionIndex)
        if HadSessions
            warning('DataIn had sessions but participant.tsv had not, skipping');
            return;
        else
            fprintf('Missing session_id column in pre-existing participant.tsv, creating\n');
            CellArrayOrig(:,3:end+1) = CellArrayOrig(:,2:end);
            CellArrayOrig{1, 2} = 'session';
            CellArrayOrig(2:end, 2) = repmat({'ASL_1'}, [size(CellArrayOrig,1)-1 1]);
        end
    end
    % 3C) Sort columns to start with participant_id & session_id
    SubjectIndex = find(cellfun(@(y) ~isempty(regexp(y,'^(participant|subject).*id*$')), lower(CellArrayOrig(1,:))));
    SessionIndex = find(cellfun(@(y) ~isempty(regexp(y,'^(session).*(id)*$')), lower(CellArrayOrig(1,:))));
    
    CellArray = CellArrayOrig(:, [SubjectIndex SessionIndex]);
    % Remove subjects & sessions from original CellArray
    CellArrayOrig = CellArrayOrig(:, [1:SubjectIndex-1 SubjectIndex+1:end]);
    SessionIndex = find(cellfun(@(y) ~isempty(regexp(y,'^(session).*(id)*$')), lower(CellArrayOrig(1,:))));
    CellArrayOrig = CellArrayOrig(:, [1:SessionIndex-1 SessionIndex+1:end]);
    CellArray = [CellArray CellArrayOrig];
    
	% If the participant ID is listed as number, then convert to char
	ParticipantIndexNumeric = find(cellfun(@(y) isnumeric(y),CellArray(1:end,1)));
	if ~isempty(ParticipantIndexNumeric)
		CellArray(ParticipantIndexNumeric,1) = cellfun(@(y) num2str(y),CellArray(ParticipantIndexNumeric,1),'UniformOutput',false);
	end

else
    % 3C) Create required columns subjects & sessions from x.SUBJECTS & x.SESSIONS

    for iSession=1:x.nSessions
        CellArray(1+iSession:x.nSessions:x.dataset.nSubjectsSessions+1, 1) = x.SUBJECTS;
    end
    CellArray(2:x.dataset.nSubjectsSessions+1, 2) = repmat(x.SESSIONS(:),[x.nSubjects, 1]);
end

% 3D) Force these names
CellArray{1,1} = 'participant_id';
CellArray{1,2} = 'session';

SubjectIndex = 1;
SessionIndex = 2;


%% -------------------------------------------------------------------------------------------
%% 4) Admin - Get column number of data
% Find columns with the same header
DataColumn = find(strcmp(CellArray(1,:), DataName));
if isempty(DataColumn)
    DataColumn = size(CellArray,2)+1;
    CellArray{1,DataColumn} = DataName;
    % now predefine/prefill with NaNs
    CellArray(2:end,DataColumn) = repmat({'n/a'}, [size(CellArray,1)-1 1]);
elseif length(DataColumn)>1
    warning('Too many columns, skipping');
    return;
end



%% -------------------------------------------------------------------------------------------
%% 5) Add data to CellArray
% Note that CellArray has headers (index +1) but DataIn has not
for iIndex=1:size(DataIn,1)
    iSubject = strcmp(CellArray(2:end,1), DataIn{iIndex,1});
    iSession = strcmp(CellArray(2:end,2), DataIn{iIndex,2});
    iRow = find(iSubject & iSession)+1; % add one for the header, assuming that DataIn has no header
    
    if isempty(iRow)
        iRow = size(CellArray,1)+1; % next row
        CellArray{iRow,1} = DataIn{iIndex,1};
        CellArray{iRow,2} = DataIn{iIndex,2};
        % Fill non-subject/non-session/non-data columns with 'n/a'
    elseif length(iRow)>1
        warning('This SubjectSession combination exists multiple times, bug? Skipping');
        continue;
    end
    
    IsEmpty = size(CellArray,1)<iRow || isempty(CellArray{iRow,DataColumn}) || strcmp(CellArray{iRow,DataColumn}, 'n/a');
    if IsEmpty || bOverwrite
        CellArray{iRow, DataColumn} = DataIn{iIndex,3};
    else
        fprintf('Data already existing, skipping\n');
    end
end


%% 6) Fill empty cells
IsEmpty = cellfun(@(y) isempty(y), CellArray);
CellArray(IsEmpty) = {'n/a'};


%% -------------------------------------------------------------------------------------------
%% 7) Sort rows on subjects
CellArray(2:end,:) = sortrows(CellArray(2:end,:), 1);

        

%% -------------------------------------------------------------------------------------------
%% 8) Write data to participants.tsv
xASL_tsvWrite(CellArray, PathTSV, 1);


end