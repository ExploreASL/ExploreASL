function [ChildJSON] = xASL_bids_InsertJSONFields(ParentJSON, ChildJSON, Fields2Skip)
%xASL_bids_InsertJSONFields InsertFields from parent JSON to child JSON
%
% FORMAT: [ChildJSON] = xASL_bids_InsertJSONFields(ParentJSON, ChildJSON[, Fields2Skip])
% 
% INPUT:
%   ParentJSON  - JSON struct or path to JSON or parms.mat (legacy) (REQUIRED)
%   ChildJSON   - JSON struct or path to JSON (legacy) (REQUIRED)
%                 This can also be a non-existing filepath, which is then
%                 created only with the provided fields in ParentJSON
%   Fields2Skip - fields from Parent JSON not to include in child JSON
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function takes all parameters from the "parent" JSON & moves them into the "child" JSON.
%               In case of co-existence of a field with different values,
%               then the value in the child JSON will prevail, per BIDS inheritance.
%
% This function runs the following steps:
%
% 1. Load JSON or parms.mat (legacy), if the inputs were paths
% 2. Insert the fields
% 3. Save a new JSON file (if ChildJSON was a path)
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:  Internally: ChildJSON = xASL_bids_InsertJSONFields(ParentJSONStruct, ChildJSONStruct, {'subject_regexp', 'ROOT'});
%           Externally: xASL_bids_InsertJSONFields('/MyStudy/ASLParameters.json', ...
%           '/MyStudy/sub-001_ses-001/perf/sub-001_ses-001_run-1_asl.json', {'subject_regexp', 'ROOT'});
% __________________________________
% Copyright 2015-2020 ExploreASL


%% 0) Admin
if nargin<3 || isempty(Fields2Skip)
    Fields2Skip = {''};
end
if nargin<2 || isempty(ChildJSON)
    error('ChildJSON input argument missing');
end
if nargin<1 || isempty(ParentJSON)
    error('ParentJSON input argument missing');
end    

PathChild = [];

%% 1) Load JSON or parms.mat
if ~isstruct(ParentJSON) && ~ischar(ParentJSON)
    error('Invalid format of ParentJSON');
elseif ~isstruct(ParentJSON) && ischar(ParentJSON)
        PathParent = ParentJSON;
    [~, Ffile, Fext] = fileparts(PathParent);
    
    if strcmp(Fext, '.json')
        % load the JSON
        ParentJSON = xASL_import_json(PathParent);
    elseif strcmp(Fext, '.mat') && strcmp(Ffile(6:end), '_parms')
        % Assume this is a parms.mat (legacy)
        mat = load(PathParent, '-mat');
        ParentJSON = mat.parms;
        % BIDS timing corrections to SI units
        if isfield(ParentJSON,'RepetitionTime')
            ParentJSON.RepetitionTime = ParentJSON.RepetitionTime/1000;
        end
        if isfield(ParentJSON,'EchoTime')
            ParentJSON.EchoTime = ParentJSON.EchoTime/1000;
        end
    else
        error('Invalid format of ParentJSON path');
    end
        
else % do nothing, assume that the struct contains JSON fields
end
if ~isstruct(ChildJSON) && ~ischar(ChildJSON)
    error('Invalid format of ChildJSON');
elseif ~isstruct(ChildJSON) && ischar(ChildJSON)
    PathChild = ChildJSON;
    % load the JSON, or initiate one
    if exist(PathChild, 'file')
        ChildJSON = xASL_import_json(PathChild);
    else
        ChildJSON = struct;
    end
else % do nothing, assume that the struct contains JSON fields
end    

%% 2) Insert fields
FieldsAre = fields(ParentJSON);

for iField=1:length(FieldsAre)
    SkipField = max(cellfun(@(y) strcmp(FieldsAre{iField},y), Fields2Skip));
    if ~SkipField
        FieldValue = ParentJSON.(FieldsAre{iField});
        if ischar(FieldValue) || isnumeric(FieldValue) || islogical(FieldValue)

            if isfield(ChildJSON,FieldsAre{iField})
                % Skip this field: per inheritance principle, daughters
                % fields have preference
                % Also skip this field if it is an environment parameter
            else
                ChildJSON.(FieldsAre{iField}) = ParentJSON.(FieldsAre{iField});
            end
        else % assume we have subfields
            try
                Subfields = fields(ParentJSON.(FieldsAre{iField}));
            catch ME
                warning(['Something went wrong with field ' FieldsAre{iField}]);
                fprintf('Is this field not char/numeric and it doesnt have subfields?\n');
                fprintf('%s\n', ME.message);
                continue; % with next field
            end
            % do the same for the subfields
            if ~isfield(ChildJSON,FieldsAre{iField})
                ChildJSON.(FieldsAre{iField}) = struct;
            end

            ChildJSON.(FieldsAre{iField}) = xASL_bids_InsertJSONFields(ParentJSON.(FieldsAre{iField}), ChildJSON.(FieldsAre{iField}));
        end
    end
end

%% 3) Save new JSON file
if ~isempty(PathChild)
    % any pre-existing JSON file will be overwritten
    spm_jsonwrite(PathChild, ChildJSON);
end


end