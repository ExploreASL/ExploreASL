function [x] = xASL_import_json(DataParFile)
% xASL_import_json This function reads in a DATA_PAR file and creates the x structure.
%
% FORMAT:   [x] = xASL_import_json(DataParFile)
%
% INPUT:
%   DataParFile     - Filename of the DATA_PAR file
%
% OUTPUT:
%   x               - x structure
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function reads in a DATA_PAR file and creates the x
%               structure. The name of the DATA_PAR file is given as a string or
%               character array. The output is the x structure.
%
%               If the DATA_PAR file is the dataset_description.json file of the BIDS
%               standard, the x structure is created according to BIDS.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:      xASL_import_json('DataParFile.json')
%
% EXAMPLE 2:    JSON FILE
%
% {
% 	"x": [{
% 			"name": 				"ExampleDataSet",
% 			"subject_regexp": 		"^Sub-\\d{3}$",
% 			"M0": 					"separate_scan",
% 			"BackgroundSuppressionNumberPulses":	"2",
% 			"readout_dim": 			"2D",
% 			"QUALITY": 				"0",
% 			"Vendor": 				"Philips",
% 			"LabelingType": 		"CASL",
% 			"qnt_init_PLD": 		"1525",
% 			"qnt_labdur": 			"1650",
% 			"qnt_PLDslicereadout": 	"43.7647"
% 		}]
% }
%
% Don't forget to escape the backslashes!
%
% OR:
%
% {
% 	"Name":         "...",
% 	"BIDSVersion":  "1.0.2",
% 	"License":      "...",
% 	"Authors":      ["..."]
% }
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2020 ExploreASL

%% Input Check
if ~exist(DataParFile, 'file')
    error('DataParFile does not exist...');
end

% Input has to be a character array
DataParFile = char(DataParFile);

%% Decode JSON file
jsonData = spm_jsonread(DataParFile);

%% Check if x field exists
if isfield(jsonData,'x')
    %% DATA_PAR file identified
    
    % Create x structure
    x = jsonData.x;
else
    x = jsonData;
end
    
% Check field names -> LEGACY
sFields = fieldnames(x);
for n=1:size(sFields,1)
    if strcmp(sFields{n},'group') % Convert group fields to correct Matlab cell arrays
        % Generate new Matlab cell array
        x.group{str2num(strrep(sFields{n},'group',''))} = x.(sFields{n});
        % Remove old field
        x = rmfield(x,sFields{n});
    elseif strcmp(sFields{n},'LOAD') % Handle load commands
        loadPaths = fieldnames(x.(sFields{n}));
        for m=1:size(loadPaths,1)
            % Don't forget to use \\ instead of \ in your paths
            if exist(x.(sFields{n}).(loadPaths{m}),'file')
                load(fullfile(x.(sFields{n}).(loadPaths{m})));
            end
        end
    end
end

%% Convert strings containing numbers to number
x = ConvertNumericalFields(x);
if isfield(x,'Q')
	x.Q = ConvertNumericalFields(x.Q);
end

%% NEW BIDS VERSION, FUTURE USE
%% elseif isfield(jsonData,'BIDSVersion')
%     %% dataset_description file identified
%     
%     % Create x structure
%     x = jsonData;
%     
%     % Get folder of dataset_description file
%     [rootBIDS,~,~] = fileparts(DataParFile);
%     if exist(rootBIDS,'dir')
%         % Check if participants file exists
%         if exist(fullfile(rootBIDS,'participants.tsv'),'file')
%             participants = tdfread(fullfile(rootBIDS, 'participants.tsv'));
%             % Check if the participants id field exists
%             if isfield(participants,'participant_id')
%                 % Get number of participants
%                 numOfP = size(participants.participant_id,1);
%                 partFieldNames = fieldnames(participants);
%                 % Convert struct to struct array, so that you can access a
%                 % participant by using: x.participant(n)
%                 for n=1:numOfP
%                     % Go trough all fields
%                     for m=1:length(partFieldNames)
%                         % Assign field to new struct
%                         x.participant(n).(partFieldNames{m})=participants.(partFieldNames{m})(n,:);
%                     end
%                 end
%             end
%         end
%         % Get task list
%         x.taskList = dir(fullfile(rootBIDS, 'task*.json'));
%         
%         % Check participants
%         if isfield(x,'participant')
%             % Get number of participants
%             numOfP = size(participants.participant_id,1);
%             % Iterate over participants
%             for n=1:numOfP
%                 % Get current participant id
%                 curID = x.participant(n).participant_id;
%                 % Check if corresponding folder exists
%                 if exist(fullfile(rootBIDS,string(curID)),'dir')
%                     % Subject file
%                     subFile = dir(fullfile(rootBIDS,string(curID),'*.tsv'));
%                     % Get session tsv file
%                     if exist(fullfile(subFile.folder,subFile.name),'file')
%                         % Display current participant
%                         disp('----------------------------------------------------------------------------------')
%                         xFieldNames = fieldnames(x.participant);
%                         for r=1:length(xFieldNames)
%                             if strcmp(xFieldNames{r},'participant_id')
%                                 fprintf('Participant:\t\t\t\t\t%s\n',string(x.participant(n).(xFieldNames{r})));
%                             elseif ~contains(xFieldNames{r},'session')
%                                 fprintf('%s:\t\t\t\t\t\t\t%s\n',xFieldNames{r},string(x.participant(n).(xFieldNames{r})));
%                             end
%                         end
%                         % Load the current session information
%                         tmpSubject = tdfread(fullfile(subFile.folder,subFile.name));
%                         fieldNames = fieldnames(tmpSubject);
%                         % Create Session Sub Structure
%                         numOfSessions = numel(fieldnames(tmpSubject));
%                         x.participant(n).session(1:numOfSessions)=struct;
%                         % Add Fields to Sub Structure
%                         for m=1:numOfSessions
%                             for p=1:length(fieldNames)
%                                 x.participant(n).session(m).(fieldNames{p})=tmpSubject.(fieldNames{p})(m,:);
%                             end
%                         end
%                         % Add Session Fields
%                         for m=1:size(x.participant(n).session,2)
%                             % Add path
%                             x.participant(n).session(m).path = fullfile(subFile.folder,x.participant(n).session(m).session_id);
%                             % Load the specific session information
%                             tableDir = dir(fullfile(x.participant(n).session(m).path,'*.tsv'));
%                             specSession = tdfread(fullfile(tableDir.folder,tableDir.name));
%                             specFieldNames = fieldnames(specSession);
%                             % Create Session Sub Structure
%                             numOfFields = numel(fieldnames(specSession));
%                             % Add Fields to Sub Structure
%                             for p=1:numOfFields
%                                 x.participant(n).session(m).(specFieldNames{p})=specSession.(specFieldNames{p});
%                             end
%                             % Display Session
%                             fprintf('Session:\t\t\t\t\t\t%s\n',x.participant(n).session(m).session_id);
%                             if isfield(x.participant(n).session(m),'filename')
%                                 for rx=1:size(x.participant(n).session(m).filename,1)
%                                     fprintf('File:\t\t\t\t\t\t\t%s\n',x.participant(n).session(m).filename(rx,:));
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end      
%         
%         % End
%         disp('----------------------------------------------------------------------------------')
%                 
%         % Work in progress ...
%         disp('Work in progress...')
%         pause;
%     end
%     
%     
% else
%     disp('Neither the old nor the new standard was used...');
% end


end

function [StructOut] = ConvertNumericalFields(StructIn)
% Also remove invalid fields    

FieldsAre = fields(StructIn);
for iField=1:length(FieldsAre)
    if length(FieldsAre{iField})>63
        warning('Invalid field name, removing from struct');
    else
        TempField = xASL_str2num(StructIn.(FieldsAre{iField}));
        if isnumeric(TempField) && size(TempField,1)>size(TempField,2)
            % enforce a horizontal vector if numerical
            TempField = TempField';
        end
        
        if min(isnumeric(TempField)) && min(~isnan(TempField))
            % if there is a numerical vector inside the string
%                 if min(isfinite(TempField))
                StructOut.(FieldsAre{iField}) = TempField;
%                 end
        else
            % keep the string
            StructOut.(FieldsAre{iField}) = StructIn.(FieldsAre{iField});
        end
    end
end


end
