function xASL_adm_BIDSJsonChecker(x,ImageType,OutputTable)
% xASL_adm_BIDSJsonChecker(x,ImageType) checks every JSON of ImageType for all subjects for differences in BIDS parameters
%
% INPUT:
%   x                            - struct containing statistical pipeline environment parameters (REQUIRED)
%   ImageType                    - string detailing JSON of which filetype is checked, default is ASL4D, other options are M0, T1 and FLAIR
%   OutputTable                  - boolean to print a table showing JSON parameters for each subject
%
% OUTPUT: n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function checks every JSON of ImageType for all subjects for differences in BIDS parameters
%
%              1. Adminstration
%                   a. Check input
%                   b. Create paths
%                   c. Create subject list
%              2. Read and save JSON parameters for each subject
%              3. Check JSON parameters for uniqueness
%              4. Save table with JSON parameter details
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_adm_BIDSJsonChecker(x,'ASL4D');
% __________________________________
% Copyright 2022 ExploreASL



%% 1. Administration
% a. Check input
if nargin <2 || isempty(ImageType)
    ImageType = 'ASL4D'; % automatically assume ASL4D json needs to be check, other options are ASL4D M0 T1 FLAIR
end

if nargin <3 || isempty(OutputTable)
    OutputTable = 0; % set writing table boolean to 0 if not provided
end

JSON = [ImageType '.json'];  % type of JSON that will be checked

% b. Create paths
TablePath = fullfile(x.dir.xASLDerivatives,'JSONcomparisonTable.tsv'); % table output

% c. Create subject list
SubjectList = xASL_adm_GetFileList(x.dir.xASLDerivatives,x.dataset.subjectRegexp,'List',[],true); % all subjects
nSubjects = numel(x.dataset.TotalSubjects);


%% 2. Read and save JSON parameters for each subject
% For loop which performs per subject creation of session list (if required), JSON path creation
% parsing of JSON's, creating of comparison table

for iSubject = 1:nSubjects
    
    % select subject and subject folder
    Subject{iSubject,1} = SubjectList{iSubject};
    SubjectFolder = fullfile(x.dir.xASLDerivatives, Subject{iSubject});
    
    % select session and session folder per subject
    Session_list = xASL_adm_GetFileList(SubjectFolder,'ASL_*','List',[],true);
    nSessions = numel(Session_list); % amount of sessions
    
    % check JSONs per subject and session
    for iSession = 1:nSessions
        if strcmp(ImageType,'ASL4D') || strcmp(ImageType,'M0')% ASL or M0 JSONs
            Session{iSession} = Session_list{iSession}; % select session
            JSONpath{iSession} = fullfile(SubjectFolder,Session{iSession},JSON); % JSON path
        else % structural JSONs
            JSONpath{iSession} = [SubjectFolder '/' JSON]; % JSON path
        end
        
        %%% read JSON
        JSONfile = spm_jsonread(char(JSONpath{iSession}));
        
        % JSON parameter removal or conversion
        if isfield(JSONfile,'GeneratedBy') % remove non-required text
            JSONfile.GeneratedBy = [];
        end
        if isfield(JSONfile,'ImageType') % convert array to char withing JSON field
            JSONfile.ImageType = char(JSONfile.ImageType);
        end
        
        SizeJSON = length(fieldnames(JSONfile)); % Nnumber of JSON parameters
        
        % create cell per subject to check for differences in JSON Fields
        JSONCell(iSubject,:) = reshape(struct2cell(JSONfile),[1 SizeJSON]);
        for iCell = 1:size(JSONCell,2)
            try
                if isnumeric(JSONCell{iCell})
                    iCellString = xASL_num2str(JSONCell{iCell});
                else
                    iCellString = char(JSONCell{iCell});
                end
                if iCell == 1
                    iSubjectJSONCellString = iCellString;
                else
                    iSubjectJSONCellString = [iSubjectJSONCellString iCellString];
                end
            catch
                warning('Could not add JSON parameter to string for checking')
            end
        end
        
        try
            JSONCellString(iSubject,:) = iSubjectJSONCellString; % add complete JSON parameter string to total
        catch
            warning('Subjects differ in ASL-BIDS JSON parameters'); % provide warning of difference between subjects
        end
        
        % put results in JSONCell containing subject and session names
        JSONCellFinal(iSubject,1) = Subject(iSubject,1);
        JSONCellFinal(iSubject,2) = Session(iSession,1);
        JSONCellFinal(iSubject,3:SizeJSON+2) = JSONCell(iSubject,:);
    end
    
end


%% 3. check all JSON fields for uniqueness for all subjects and sessions
for i = 1:size(JSONCell,2)
    if i > 2 % skip subject and session
        try UniqueResult = unique(string(JSONCell(:,i)),'rows'); % check for more than one unique parameter value
            if size(UniqueResult,1) > 1
                warning(['JSON parameter ' JSONtableFinale.Properties.VariableNames{i+2} ' is different for subjects']) % add 2 for difference in table and cell (Subject name and session name)
            end
        end
    end
end


%% 4. write Cell to .tsv
if OutputTable == 1
    xASL_tsvWrite(JSONCellFinal,TablePath,1,0) 
end
end

