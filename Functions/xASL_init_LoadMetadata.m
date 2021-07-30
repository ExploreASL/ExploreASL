function [x] = xASL_init_LoadMetadata(x)
%LoadMetadata Load all study metadata
%
% FORMAT: [x] = xASL_init_LoadMetadata(x)
%
% INPUT:
%   x                   - struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
% INPUT FILES:
% participants.tsv      - per BIDS, file in root analysis folder containing
%                         the same as the individual .mat files explained
%                         below
% [VariableName].mat    - ExploreASL legacy, each .mat file contains a
%                         different metadata variable. Any string or
%                         numerical value is possible. For missing data,
%                         use 'n/a' (per BIDS).
%                         Subject and session/run names should be unique.
% OUTPUT:
%   x   - struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%   x.S.SetsName        - name of loaded variable (cell array)
%   x.S.SetsID          - metadata content (numerical, double. 0 stands for missing data)
%   x.S.Sets1_2Sample   - Contains the sample datatype (vector), options are:
%                         1) ordinal data, treat as paired. Allows
%                         comparing between groups, e.g. session/run,
%                         timepoint/visit;
%                         2) ordinal data, treated as unpaired. Allows
%                         statistically comparing between groups or accounting for group
%                         differences in image processing, e.g.
%                         Site/Sequence;
%                         3) continuous data, e.g. age.
%                         All data is assumed to be continuous, unless
%                         few unique options (set by MaxNumOrdinalOptions below)
%                         or explicitly stated elsewhere.
%   x.S.SetsOptions     - Name of options of ordinal data (cell array of strings), e.g. "Male"
%                         "Female" for sex, e.g. "AmsterdamUMC"
%                         "UniversityCollegeLondon" for Site.
%                         For continuous data, this is ignored. Only one
%                         option will be specified, with string equal to
%                         x.S.SetsName. Note also that ordinal data with
%                         options "0" and "1", will be kept like this. 
%
%                         E.g. sex with input options "0" and "1" or "male"
%                         "female" will both be coded in x.S.SetsID as 1
%                         and 2 (0 for missing data). 
%                         x.S.SetsOptions will be {"0" "1"} or {"male"
%                         "female"} respectively. 
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function loads all metadata used in the study, either statistical
% covariates (age, MMSE) or groups to compare between (site, sequence,
% cohort), or parameters to be used in quantification/image processing
%
% These parameters should be provided in .mat files in the root analysis
% folder. Each .mat file should contain a single type of metadata, and 
% the filename should equal the variable name.
% Metadata content should be a cell array
% with subjects as first column and metadata as last column.
% Sessions (runs) can be included as second column.
%
% Metadata can be in any string or numerical format.
%
% participants.tsv is now added per BIDS. It looks for metadata in participants.tsv first,
% before going through the mat-files
%
% This function iterates through the following steps for each variable:
%
% 1. Admin (what nOptions do we call ordinal, convert subject numeric to string, remove white spaces from data)
% 2. Get unique list of data options & check for missing data
% 3. Deal with data format (correct NaNs, deal with numeric vs strings)
% 4. Distinguish continous data (e.g. age) or ordinal data (groups to compare, e.g. cohort)
% 5. Check if data is complete for all subjects
% 6. Include complete data in x.S.SETS
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: x = xASL_init_LoadMetadata(x);
% __________________________________
% Copyright 2015-2020 ExploreASL




    %% ------------------------------------------------------------------------------------------------------------
    %% Admin
    if x.nSubjects==0
        fprintf('%s\n','No variables loaded because no subjects found');
        return;
    end

    
    
    x = xASL_bids_LoadParticipantTSV(x);
    
    if isfield(x,'S') && isfield(x.S, 'SetsName')
        SetsInParticipantTSV = x.S.SetsName;
    else
        SetsInParticipantTSV = '';
        fprintf('No participants.tsv file detected\n');
    end
    
    MatFileList = xASL_adm_GetFileList(x.D.ROOT, '^(?!(x|xASL|dcm2niiCatchedErrors)).*\.mat$', 'FPList', [0 Inf]);
    
    x = xASL_init_LoadMat(x, MatFileList); % legacy
    
    %% 
    %% -----------------------------------------------------------------------------------------------
    %% Add stats in participants.tsv
    bCreateParticipantsTsv = 0;
    
    if ~isempty(MatFileList)
        warning('Legacy .mat-files detected in analysis rootfolder');
        if ~isfield(x, 'S') || ~isfield(x.S, 'SetsName')
            fprintf('But not added as sets\n');
            fprintf('Were these .mat files without participant information?\n');
            fprintf('Any participant information inside a .mat-file should be in a variable with identical name as the .mat filename\n');
        elseif isempty(SetsInParticipantTSV)
            fprintf('Creating participants.tsv file from them\n');
            fprintf('Consider deleting the mat-files and keep the participants.tsv file only\n');
            bCreateParticipantsTsv = 1;
        elseif isequal(SetsInParticipantTSV, x.S.SetsName)
            warning('Both participants.tsv & legacy .mat-files existed, with the same parameters');
            fprintf('Please merge these manually, and/or delete any old .mat-files\n');
        else
            warning('Both participants.tsv & legacy .mat-files existed, trying to merge them');
            fprintf('Consider deleting the mat-files and keep the participants.tsv file only\n');
            bCreateParticipantsTsv = 1;
        end
    end
        
    if bCreateParticipantsTsv
        for iSession=1:x.dataset.nSessions
            VarDataOri([iSession:x.dataset.nSessions:x.dataset.nSubjectsSessions-x.dataset.nSessions+iSession], 1) = x.SUBJECTS(:);
        end
        VarDataOri(:,2) = repmat(x.SESSIONS(:), [x.nSubjects 1]);

        for iSet=1:length(x.S.SetsName)
            % skip creating metadata columns for some parameters
            % (session/subjectnlist)
            if isempty(regexpi(x.S.SetsName{iSet},'(session|subjectnlist)'))
                % Fill VarData
                VarData = VarDataOri;

                for iSubjectSession=1:size(x.S.SetsID,1)
                    if x.S.Sets1_2Sample(iSet)==3
                        % x.S.SetsID contains continuous numbers
                        % x.S.SetsOptions == x.S.SetsName
                        VarData{iSubjectSession,3} = x.S.SetsID(iSubjectSession,iSet);
                    else % x.S.SetsID contains category number
                        % x.S.SetsOptions has the name of the category
                        if isfinite(x.S.SetsID(iSubjectSession,iSet))
                            VarData{iSubjectSession,3} = x.S.SetsOptions{iSet}{x.S.SetsID(iSubjectSession,iSet)};
                        else
                            VarData{iSubjectSession,3} = 'n/a';
                        end
                    end
                end

                xASL_bids_Add2ParticipantsTSV(VarData, x.S.SetsName{iSet}, x, 1); % overwrite
            end
        end
    end
        
    %% ------------------------------------------------------------------------------------------------------------


    
end



function [x] = xASL_bids_LoadParticipantTSV(x)
%xASL_bids_LoadParticipantTSV 

    %% Load participants.tsv
    PathTSV = fullfile(x.D.ROOT, 'participants.tsv');
    bParticipantsTSV = false;

    % Backwards compatibility: rename to lowercase
    fileListParticipantsTSVold = xASL_adm_GetFileList(x.D.ROOT,'^Participants.tsv$',false);
    if ~isempty(fileListParticipantsTSVold) % We added the _temp copy step so that the code works on case insensitive systems like windows as well. Please don't remove that step for backwards compatibility (at least not until release 2.0.0).
        pathParticipantsTSVold = fullfile(x.D.ROOT, 'Participants.tsv');
        pathParticipantsTSVoldTemp = fullfile(x.D.ROOT, 'Participants_temp.tsv');
        xASL_Move(pathParticipantsTSVold,pathParticipantsTSVoldTemp);
        xASL_Move(pathParticipantsTSVoldTemp,PathTSV);
    end

    % Read in cell array
    if exist(PathTSV, 'file')
        CellArray = xASL_tsvRead(PathTSV);
        bParticipantsTSV = true;
    end

    if bParticipantsTSV
        fprintf('participants.tsv (BIDS) detected, loading...\n');
        
        % Check if the TSV has empty cells, because this can cause crashes & is illegal (per BIDS)
        HasEmptyCells = false;
        for iX=1:size(CellArray, 1)
            for iY=1:size(CellArray, 2)
                if isempty(CellArray{iX, iY})
                    HasEmptyCells = true;
                end
            end
        end
        if HasEmptyCells
            warning('participants.tsv contains empty cells, missing data should be filled by n/a (per BIDS)');
        end           
            
        
        %% First isolate fixed variables (subject, session/run, site)
        SubjectIndex = find(cellfun(@(y) ~isempty(regexp(y,'^(participant|subject).*id*$')), lower(CellArray(1,:))));
        SessionIndex = find(cellfun(@(y) ~isempty(regexp(y,'^(session).*(id)*$')), lower(CellArray(1,:))));
        SiteIndex = find(cellfun(@(y) ~isempty(regexp(y,'^(site).*(id)*$')), lower(CellArray(1,:))));
        
        if isempty(SubjectIndex)
            warning('Couldnt find subject index, skipping this participant.tsv');
            return;
        end
        
        HasSession = ~isempty(SessionIndex);
        HasSite = ~isempty(SiteIndex);
        
        %% Separate data from subjectIDs
        SubjectColumn = CellArray(2:end, SubjectIndex);
        CellArray = CellArray(:,[1:SubjectIndex-1 SubjectIndex+1:end]);
        
        if HasSession
            % Separate data from sessionIDs
            % First obtain new SessionIndex
            SessionIndex = find(cellfun(@(y) ~isempty(regexp(y,'^(session).*(id)*$')), CellArray(1,:)));
            SessionColumn = CellArray(2:end, SessionIndex);
            CellArray = CellArray(:, [1:SessionIndex-1 SessionIndex+1:end]);
        end

        
        %% Sort site column as first column
        if HasSite
            SiteIndex = find(cellfun(@(y) ~isempty(regexp(y,'^(site).*(id)*$')), lower(CellArray(1,:))));
            CellArray = CellArray(:, [SiteIndex 1:SiteIndex-1 SiteIndex+1:end]);
        end
        
        %% Skip if no data
        
        
        %% Loop over metadata columns/variables
        for iVar=1:size(CellArray,2)
            VarName = CellArray{1, iVar};
            VarOptions = {VarName}; % AddMetadata function will automatically detect if unique options can be defined
            VarSample = 2; % AddMetadata function will automatically detect if otherwise
            VarContent = SubjectColumn;
            if HasSession
                VarContent(:,2) = SessionColumn;
            end
            VarContent(:,end+1) = CellArray(2:end, iVar);
            x = xASL_init_AddVariable(x, VarName, VarContent, VarOptions, VarSample);
        end
    end

end


%% ------------------------------------------------------------------------------------------------------------
function [x] = xASL_init_LoadMat(x, MatFileList)
%xASL_init_LoadMat

    %% Add from mat-files. Start with site
    % Here we loop over .mat files to load as covariates/parameters to
    % be used in the pipeline (e.g. quantification parameters such as
    % hematocrit/bloodT1/T2star, or study parameters such as timepoint
    % (visit), session (run), cohort etc.
    SiteMat = fullfile(x.D.ROOT, 'Site.mat');

    if exist(SiteMat,'file')
        FList{1} = SiteMat;
        FList(2:1+length(MatFileList),1) = MatFileList;
    else
        FList = MatFileList;
    end

    for iL=1:length(FList)
        [~, FileName] = fileparts(FList{iL});

        VarName = FileName;
        VarOptions = {FileName};
        VarSample = 2; % (datatype, 0=one-sample, 1=paired, 2=two-sample, 3=continuous
        % Here, we assume that the datatype is two-sample. The script below
        % "AddStatisticalVariable" will automatically change this into a
        % continuous variable if length(unique)>5

        MatContent = load(FList{iL});
        if isfield(MatContent,FileName)
            VarContent  = MatContent.(FileName);
            x = xASL_init_AddVariable(x, VarName, VarContent, VarOptions, VarSample);
        else
            fprintf('%s\n',['mat-file ' FileName '.mat did not contain stats variable equal to variable name & was skipped']);
        end
    end


end



%% ------------------------------------------------------------------------------------------------------------
function [x] = xASL_init_AddVariable(x, VarName, VarContent, VarOptions, VarSample)
%xASL_init_AddVariable Load, process & add a single variable to memory
%
% VarName       - equal to x.S.SetsName specified in the function header
% VarContent    - equal to x.S.SetsID
% VarOptions    - equal to x.S.SetsOptions
% VarSample     - equal to x.S.Sets1_2Sample


%% ------------------------------------------------------------------------------------------------------------
%% 1) Admin
% A) Max num ordinal options
% B) If subject name are numbers, convert to strings
% C) Remove white spaces from data

% Maximum number of ordinal options. This makes sure that if there are
% more options than this number, and variate was loaded as "ordinal", it
% is set to be "continuous"

MaxNumOrdinalOptions = 15;
IfProcess = 1;

if  size(VarContent,2)<2
    fprintf('%s\n',['Variable ' VarName ' skipped, did not contain data or only subject names']);
    IfProcess= 0; % skip this variable, doesn't contain information
end    

% Check if already exists, then remove previous set
if isfield(x.S,'SetsName') && IfProcess
    iSet = find(strcmp(x.S.SetsName,VarName));
    if ~isempty(iSet)
        % remove existing set, to add new one
        AllSets = 1:length(x.S.SetsName);
        AllSets = AllSets([1:iSet-1 iSet+1:end]);
        x.S.SetsName = x.S.SetsName(AllSets);
        x.S.SetsID = x.S.SetsID(:,AllSets);
        x.S.SetsOptions = x.S.SetsOptions(AllSets);
        x.S.Sets1_2Sample = x.S.Sets1_2Sample(AllSets);
    end
end

if ~IfProcess
    % if variable is invalid to process (e.g. already exist)
    return;
end

% A) Find SetIndex
if ~isfield(x.S,'SetsName')
    SetIndex = 1;
else
    SetIndex = length(x.S.SetsName)+1;
end

if ~iscell(VarContent)
    warning('Metadata was not defined as cell array');
    fprintf('%s\n', ['Skipping variable ' VarName]);
    return;
end

% B) If SubjectNames are numbers, convert to strings
for iNum=1:size(VarContent,1)
    if isnumeric(VarContent{iNum,1})
        VarContent{iNum,1}  = num2str(VarContent{iNum,1});
    end
end


% C) remove whitespaces
VarContent(:,1) = strtrim(VarContent(:,1));

% Determine whether group contains subject-wise or 
% SubjectSession-wise values
ContainsSessions = xASL_adm_FindStrIndex( VarContent,x.SESSIONS{1});
if ~isempty(ContainsSessions) && min(min(ContainsSessions~=0))
    ContainsSessions = 1;
    SessionColumn = 2;
    DataColumn = 3;
else
    ContainsSessions = 0;
    DataColumn = 2;
end




%% ------------------------------------------------------------------------------------------------------------
%% 2) Get unique list of data options & check for empty lines
% Either for subjects only, or when including sessions
nextN=1;
fprintf('%s\n',['Loading ' VarName ':  ']);
for iS=1:x.nSubjects
    xASL_TrackProgress(iS,x.nSubjects);

    HasEmpty = 0;
    IsFound = 0;
    SubjectIndex = find(strcmp(VarContent(:,1), x.SUBJECTS{iS}));

    if isempty(SubjectIndex)
        % retry
        for iiG=1:size(VarContent,1)
            if strcmp(VarContent{iiG,1}, x.SUBJECTS{iS})
                SubjectIndex = iiG;
            end
        end
    end
    
    
    for iSess=1:x.dataset.nSessions
        if DataColumn==3 && length(SubjectIndex)>1
            % if there are sessions, account for them
            iG = find(strcmp(VarContent(SubjectIndex,2),x.SESSIONS{iSess}));
            iG = SubjectIndex(iG);
        else
            iG = SubjectIndex;
        end

        if ~isempty(iG)
            IsFound = 1;
            iG = iG(1);
            % if subjects or sessions are mentioned multiple times
            % in the DataFile, use only the first instance                    

            if isempty(VarContent{iG,DataColumn}) % check first if it is not empty
                fprintf('%s\n',['Subject #' num2str(iG) ', ' VarContent{iG,1} ' is empty for variable ' VarName]);
                HasEmpty = 1;
            else
                TempDataList{nextN} = VarContent{iG,DataColumn};

                if islogical(TempDataList{nextN}) % correct for logical
                    if TempDataList{nextN}==0
                        TempDataList{nextN} = 0;
                    else
                        TempDataList{nextN} = 1;
                    end
                end
                nextN = nextN+1;
            end
        end
        if ~IsFound; HasEmpty = 1; end;
    end
end

fprintf('\n');

if ~exist('TempDataList','var') || isempty(TempDataList)
    fprintf('%s\n',['Variable ' VarName ' not included because all subjects were missing']);
    return;
end



% If data is a combination of strings & numbers, convert all to strings
CountString = 0;
CountNumber = 0;

for iN=1:length(TempDataList)
    if isstr(TempDataList{iN})
            CountString = CountString+1;
    elseif  isnumeric(TempDataList{iN})
            CountNumber = CountNumber+1;
    end
end

if CountString>0 && CountNumber>0 % if both strings & number exists
    % convert all to strings
    for iN=1:length(TempDataList)
        TempDataList{iN} = num2str(TempDataList{iN});
    end
end

try 
    DataOptionList = unique(TempDataList);
catch
    % this workaround is here because "unique" doesnt accept
    % cells with numbers, only numbers or cells with strings
    for iD=1:length(TempDataList)
        NumericDataList(iD) = TempDataList{iD};
    end
    DataOptionList = unique(NumericDataList);
end

DataOptionList = sort(DataOptionList);




%% ------------------------------------------------------------------------------------------------------------
%% 3) Deal with data format

% A) Correct DataOptionList for multiple NaNs, if data is numeric
if isnumeric(DataOptionList)
    NextN = 1;
    NaNindex = 1;
    for iD=1:length(DataOptionList)
        if isnan(DataOptionList(iD))
            NaNindex(NextN) = iD;
            NextN = NextN+1;
        end
    end

    if sum(NaNindex>1)
        DataOptionList = DataOptionList(1:min(NaNindex));
    end        
end



% % B) If data are numbers hidden in strings, convert to numbers
% %    Else, convert strings to ordinal numbers
% for iNum=1:size(VarContent,1)
%     if  islogical(VarContent{iNum,DataColumn})
%         if  VarContent{iNum,DataColumn}==0
%             VarContent{iNum,DataColumn} = 0;
%         else
%             VarContent{iNum,DataColumn} = 1;
%         end
%     elseif ~isnumeric(VarContent{iNum,DataColumn})
%         try
%             VarDataIs = str2num(VarContent{iNum, DataColumn});
%         catch ME
%             fprintf('%s\n', ME.message);
%             error(['Covariant data of ' VarName ' were not numbers?']);
%             
%         end
% 
%         if  isempty(VarDataIs)
%             % Assume that data are strings, so take unique strings & assign numbers to them
%             VarDataIs = find(strcmp(DataOptionList,VarContent{iNum,DataColumn}));
%          end
% 
%         try
%             VarContent{iNum,DataColumn} = VarDataIs;
%         catch ME
%             fprintf('%s\n', ME.message);
%             error(['Covariant data of ' VarName ' were not numbers?']);
%         end                
%     end
% end






%% ------------------------------------------------------------------------------------------------------------
%% 4) Distinguish continous data (e.g. age) or ordinal data (groups to compare, e.g. cohort)
%% 4A) Force continuous data for known parameters
ListContinuousData = {'AcquisitionTime' 'age' 'MeanMotion' 'GM_vol' 'WM_vol' 'CSF_vol' 'GMWM_ICVRatio' 'GM_ICVRatio' 'WMH_count' 'WMH_vol'...
    'CBF_spatial_CoV' 'CBF_spatial_CoV_norm' 'PseudoCBF_spatial_CoV' 'Yrs_AAO' 'Hematocrit' 'SliceReadoutTime' 'AcquisitionTime' 'Inititial_PLD'...
    'LabelingEfficiency' 'qnt_ATT' 'qnt_T1a' 'qnt_lab_eff'...
    'Blood_T1art' 'GM','deep_WM','WholeBrain','L-ICA','R-ICA','POS','Caudate','Cerebellum','Frontal','Insula','Occipital','Parietal','Putamen'...
    'Temporal','Thalamus','ACA_1','ACA_2','ACA_3','MCA_1','MCA_2','MCA_3','PCA_1','PCA_2','PCA_3','ACA_bilat','MCA_bilat','PCA_bilat','WMH','NAWM'};

% this is a list of known variables that will never be ordinal, always continuous
% Check if current VarName is one of these

DataIsContinuous = ~isempty(find(strcmp(lower(ListContinuousData),lower(VarName))));
DataIsContinuous = length(DataOptionList)>MaxNumOrdinalOptions || DataIsContinuous;

if DataIsContinuous
    VarSample = 3; % continuous
    x.S.Sets1_2Sample(SetIndex) = 3;
end


%% 4B) If data are considered ordinal, convert data to indices
if length(DataOptionList)<MaxNumOrdinalOptions && ~isempty(DataOptionList) && ~DataIsContinuous
    %% Specify variable options, only if the variate is ordinal
    if isnumeric(DataOptionList)
        for iD=1:length(DataOptionList)
            VarOptions{iD} = num2str(DataOptionList(iD));
        end
    else
        VarOptions = DataOptionList;
    end

    VarSample = 2; % independent samples
    
    % If data are categorical numbers, but still not ordinal numbers, convert to ordinal numbers
    % This is required because only integer numbers>0 can be indices
    % This will also avoid zeros, since 0's are no valuable index
    if ~isnumeric(DataOptionList) || ~min(DataOptionList==1:length(DataOptionList)) % if numbers are not continuous integers (== ordinal sequence)
        if isnumeric(DataOptionList)
            % we convert the option list to cell
            CurrentOptionList = num2cell(DataOptionList);
        else
            % we convert any string values to numeric if they are numerical
            CurrentOptionList = xASL_str2num(DataOptionList, 1, 0);
        end        
        
        for iNum=1:size(VarContent,1)
            if isempty(VarContent{iNum,DataColumn})
                % do nothing, this will stay empty
            else % find the corresponding DataOption, to provide this as index
                CurrentVarContent = VarContent{iNum,DataColumn};
                if ~isnumeric(CurrentVarContent) && ~isnan(xASL_str2num(CurrentVarContent))
                    % convert to numerical if appropriate
                    CurrentVarContent = xASL_str2num(CurrentVarContent);
                end
                
                if isnumeric(CurrentVarContent)
                    NumericalIndices = find(cellfun(@(y) isnumeric(y), CurrentOptionList));
                    if isnan(CurrentVarContent)
                        IndexIs = find(cellfun(@(y) isnan(y), CurrentOptionList(NumericalIndices)));
                    else
                        IndexIs = find(cellfun(@(y) CurrentVarContent==y, CurrentOptionList(NumericalIndices)));
                    end
                    CurrentVarContent = NumericalIndices(IndexIs);
                elseif ischar(CurrentVarContent)
                    StringIndices = find(cellfun(@(y) ischar(y), CurrentOptionList));
                    IndexIs = find(cellfun(@(y) strcmp(CurrentVarContent, y), CurrentOptionList(StringIndices)));
                    CurrentVarContent = StringIndices(IndexIs);
                else
                    CurrentVarContent = [];
                end
                VarContent{iNum,DataColumn} = CurrentVarContent;
            end
        end
    end
end




%% ------------------------------------------------------------------------------------------------------------
%% 5) Check if data is complete for all subjects
CountAbsent = 0;

for iSubject=1:x.nSubjects
    for iSess=1:x.dataset.nSessions

        % ID (which name, group etc), all for identification
        iSubjSess = (iSubject-1)* x.dataset.nSessions +iSess;

        % For each group, create an additional set
        % NB: at least all subjects/sessions that are used should be 
        % present, otherwise this will give an error!
        % User is responsible for imputing/interpolating missing values

        % 1st column should be subject-id column
        % 2nd column should be parameter values
        % or 2nd column is session-id (e.g. 'ASL_1', then
        % 3rd column contains parameter values


        NewVariable(iSubjSess,1) = NaN; % pre-filling for missing values
        % SubjectIndex
        SubjectIndex = xASL_adm_FindStrIndex(VarContent(:,1), x.SUBJECTS{iSubject});

        if max(isempty(SubjectIndex) | SubjectIndex==0) % if this subject couldn't be found in variable data
            AbsentSubjects{CountAbsent+1,1} = x.SUBJECTS{iSubject};
            CountAbsent = CountAbsent+1;                    

        else
            if ~ContainsSessions % If there is no session column
                VC = VarContent{SubjectIndex,DataColumn}; % this debugs any erroneous double entries
                if ~isempty(VC)
                    NewVariable(iSubjSess,1) = VC(1);
                end

            else % if there is a session column

                VarContent(:,2) = strtrim(VarContent(:,2)); % remove whitespaces

                % This code assumes 1 column for x.P.SubjectID & 1 column for x.P.SessionID SessionIndex
                SessionIndex = xASL_adm_FindStrIndex( VarContent(SubjectIndex(:,1),SessionColumn), x.SESSIONS{iSess} );

                if isempty(SessionIndex) || SessionIndex==0 % no session data found for this subject/session
                    AbsentSubjects{CountAbsent+1,1} = x.SUBJECTS{iSubject};
                    CountAbsent = CountAbsent+1;
                else
                    CurrentVarContent = xASL_str2num(VarContent{ SubjectIndex(SessionIndex),DataColumn});
                    if ~isnan(CurrentVarContent) && ~isempty(CurrentVarContent)

                        nIndex = SubjectIndex(SessionIndex);
                        VC = VarContent{nIndex, DataColumn}; % this debugs any erroneous double entries
                        if ~isempty(VC)
                            NewVariable(iSubjSess,1) = VC(1);
                        end
                    else
                        AbsentSubjects{CountAbsent+1,1} = x.SUBJECTS{iSubject};
                        CountAbsent = CountAbsent+1;
                    end % if ~isnan(CurrentVarContent) && ~isempty(CurrentVarContent)
                end % if isempty(index2) | index2==0
            end % if ~ContainsSessions
        end % if  max(isempty(index1) | index1==0
    end % for iSess=1:x.dataset.nSessions
end % for iSubj=1:x.nSubjects





%% ------------------------------------------------------------------------------------------------------------
%% 6) Include complete data in x.S.SETS
%  Incomplete data is fine, when it is continuous data, to be filled by NaNs
NewVariableUnique = unique(NewVariable(~isnan(NewVariable))); % get unique values, excluding NaNs
ManyAbsent = CountAbsent>0.05*x.nSubjects || HasEmpty;
CheckLength = length(NewVariableUnique)<0.25*x.nSubjects;
CheckLength2 = size(NewVariable, 1)~=size(x.S.SetsID, 1);

if (ManyAbsent && CheckLength) || CheckLength2 % don't include this variable because of missing data 
    % Allow for 5% missing data, if data is continuous
    % (assuming continuous data has at least 25% of sample size
    % as unique values
    if exist('AbsentSubjects', 'var')
        warning('%s\n',['Variable ' VarName ': following subjects were missing & set to NaN:']);
        for iAb=1:size(AbsentSubjects,1)
            fprintf('%s',[AbsentSubjects{iAb,1} ', ']);
            if  (iAb/12)==ceil(iAb/12) % 12 subjects per line
                fprintf('\n'); % next line
            end
        end
    else
        warning('%s\n',['Variable ' VarName ': something wrong with this variable, data missing']);
    end
    fprintf('\n');
end

if size(NewVariable, 1)==size(x.S.SetsID, 1) % include variable
    x.S.SetsID(:,SetIndex) = NewVariable;

    % Add variable content, name, option & sets_1_2_sample
    x.S.SetsName{SetIndex} = VarName;
    x.S.Sets1_2Sample(SetIndex) = VarSample;         
    x.S.SetsOptions{SetIndex} = VarOptions;

%% DELETE THIS PART LATER, SEEMS COPY OF ABOVE
%     % if data seems continuous, force Sets1_2Sample==3
%     ListContinuousData = {'AcquisitionTime' 'spatial_CoV' 'ICVRatio' 'age' 'MeanMotion' 'GM_vol' 'WM_vol' 'CSF_vol' 'GMWM_ICVRatio' 'GM_ICVRatio' 'WMH_count' 'WMH_vol'...
%         'CBF_spatial_CoV' 'CBF_spatial_CoV_norm' 'PseudoCBF_spatial_CoV' 'Yrs_AAO' 'Hematocrit' 'SliceReadoutTime' 'AcquisitionTime' 'Inititial_PLD'...
%         'Blood_T1art' 'GM','deep_WM','WholeBrain','L-ICA','R-ICA','POS','Caudate','Cerebellum','Frontal','Insula','Occipital','Parietal','Putamen'...
%         'Temporal','Thalamus','ACA_1','ACA_2','ACA_3','MCA_1','MCA_2','MCA_3','PCA_1','PCA_2','PCA_3','ACA_bilat','MCA_bilat','PCA_bilat','WMH','NAWM'};
% 
%     % this is a list of known variables that will never be ordinal, always continuous
%     % These are common continuous variables generated by ExploreASL, and should be treated as continuous
%     % Check if current VarName is one of these
%     DataIsContinuous = ~isempty(find(strcmp(lower(ListContinuousData),lower(VarName))));                
%     UniqueValues = x.S.SetsID(:,SetIndex);
%     UniqueValues = unique(UniqueValues(isfinite(UniqueValues)));
%     
%     if DataIsContinuous || (length(UniqueValues)>MaxNumOrdinalOptions && x.S.Sets1_2Sample(SetIndex)~=3)
%         x.S.Sets1_2Sample(SetIndex) = 3;
%     end



%% DELETE THIS PART LATER, SEEMS COPY OF ABOVE
    % ------------------------------------------------------------------------------------------------------------
    % Add ordinal variable options, if they are not yet there
    if x.S.Sets1_2Sample(SetIndex)~=3 && length(x.S.SetsOptions{SetIndex})==1
        UniqueData = unique( x.S.SetsID(:,SetIndex) );
        for iOption=1:length(UniqueData)
            x.S.SetsOptions{SetIndex}{iOption} = num2str(UniqueData(iOption));
        end
    end
end
   

end
