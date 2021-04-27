function jsonOut = xASL_bids_JsonCheck(jsonIn,fileType)
%xASL_bids_JsonCheck Goes through the JSON structure before saving it and ensures that it fits the ASL-BIDS format
%
% FORMAT: jsonOut = xASL_bids_JsonCheck(jsonIn,fileType)
%
% INPUT:
%   jsonIn  - JSON with the input fields (REQUIRED)
%   fileType   - 'ASL' or 'M0' when the input and output is ASL-BIDS, empty for normal BIDS (REQUIRED)
%
% OUTPUT: 
%   jsonOut - ordered and checked JSON structure
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:
% It checks the existence of all BIDS fields, removes superfluous fields, checks all the conditions and orderes
% the structure on the output. It works according to the normal BIDS, or ASL-BIDS definition
%
% EXAMPLE: n/a
%
% __________________________________
% Copyright 2015-2020 ExploreASL

% Create an empty output structure and a structure with fields to delete
jsonOut = struct;
jsonRemove = struct;

bidsPar = xASL_bids_Config();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove fields not belonging to BIDS

% Remove certain empty fields
for iField = 1:length(bidsPar.listRemoveIfEmpty)
    if isfield(jsonIn,bidsPar.listRemoveIfEmpty{iField})
        if isempty(jsonIn.(bidsPar.listRemoveIfEmpty{iField}))
            jsonRemove.(bidsPar.listRemoveIfEmpty{iField}) = '';
        end
    end
end

% Remove non-BIDS fields
for iField = 1:length(bidsPar.listFieldsRemoveGeneral)
    if isfield(jsonIn,bidsPar.listFieldsRemoveGeneral{iField})
        jsonRemove.(bidsPar.listFieldsRemoveGeneral{iField}) = '';
    end
end

% Remove non-ASL-BIDS fields from ASL sequences
if strcmpi(fileType,'ASL') || strcmpi(fileType,'M0')
    for iField = 1:length(bidsPar.listFieldsRemoveASL)
        if isfield(jsonIn,bidsPar.listFieldsRemoveASL{iField})
            jsonRemove.(bidsPar.listFieldsRemoveASL{iField}) = '';
        end
    end
else % And remove certain fields only from non-ASL sequences
	for iField = 1:length(bidsPar.listFieldsRemoveNonASL)
		if isfield(jsonIn,bidsPar.listFieldsRemoveNonASL{iField})
			jsonRemove.(bidsPar.listFieldsRemoveNonASL{iField}) = '';
		end
	end
end

% Remove fields belonging to dataset_description
for iField = 1:length(bidsPar.datasetDescription.Required)
    if isfield(jsonIn,bidsPar.datasetDescription.Required{iField})
        jsonRemove.(bidsPar.datasetDescription.Required{iField}) = '';
    end
end
for iField = 1:length(bidsPar.datasetDescription.Recommended)
    if isfield(jsonIn,bidsPar.datasetDescription.Recommended{iField})
        jsonRemove.(bidsPar.datasetDescription.Recommended{iField}) = '';
    end
end
for iField = 1:length(bidsPar.datasetDescription.Optional)
    if isfield(jsonIn,bidsPar.datasetDescription.Optional{iField})
        jsonRemove.(bidsPar.datasetDescription.Optional{iField}) = '';
    end
end

% Go through all input fields and copy to output, but skip those in jsonRemove
for nameField = fieldnames(jsonIn)'
    if ~isfield(jsonRemove, nameField{1})
        jsonOut.(nameField{1}) = jsonIn.(nameField{1});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check field requirements and dependencies

% Check required ASL fields
if strcmpi(fileType,'ASL')
    strReport = '';
    for iField = 1:length(bidsPar.ASLfields.Required)
        if ~isfield(jsonOut,bidsPar.ASLfields.Required{iField})
            if isempty(strReport)
                strReport = bidsPar.ASLfields.Required{iField};
            else
                strReport = [strReport ', ' bidsPar.ASLfields.Required{iField}];
            end
        end
    end
    if ~isempty(strReport)
        fprintf('%s\n\n',['Missing required ASL fields: ' strReport]);
    end
    
    strReport = '';
    for iField = 1:length(bidsPar.ASLfields.Recommended)
        if ~isfield(jsonOut,bidsPar.ASLfields.Recommended{iField})
            if isempty(strReport)
                strReport = bidsPar.ASLfields.Recommended{iField};
            else
                strReport = [strReport ', ' bidsPar.ASLfields.Recommended{iField}];
            end
        end
    end
    if ~isempty(strReport)
        fprintf('%s\n\n',['Missing Recommended ASL fields: ' strReport]);
    end
end

% Check required M0 fields
if strcmpi(fileType,'M0')
    strReport = '';
    for iField = 1:length(bidsPar.M0fields.Required)
        if ~isfield(jsonOut,bidsPar.M0fields.Required{iField})
            if isempty(strReport)
                strReport = bidsPar.M0fields.Required{iField};
            else
                strReport = [strReport ', ' bidsPar.M0fields.Required{iField}];
            end
        end
    end
    if ~isempty(strReport)
        fprintf('%s\n\n',['Missing required M0 fields: ' strReport]);
    end
end

% Check ASL dependencies
if strcmpi(fileType,'ASL')
    for iCond = 1:length(bidsPar.ASLCondition)
        % First check if the field is present
        if isfield(jsonOut,bidsPar.ASLCondition{iCond}.field)
            
            % Checking if the condition is met, assuming no
            bCond = 0;
            if isempty(bidsPar.ASLCondition{iCond}.value)
                % Empty value means the field is only present
                bCond = 1;
            elseif ischar(bidsPar.ASLCondition{iCond}.value)
                % strings are checked with regexpi
                if regexpi(jsonOut.(bidsPar.ASLCondition{iCond}.field),bidsPar.ASLCondition{iCond}.value)
                    bCond = 1;
                end
            elseif isequal(jsonOut.(bidsPar.ASLCondition{iCond}.field),bidsPar.ASLCondition{iCond}.value)
                % Logical and numbers are checked with isequal
                bCond = 1;
            end
            
            % Conditions are met, now check for dependencies
            if bCond
                % Check the required filled fields
                strReportFilled = '';
                for iField = 1:length(bidsPar.ASLCondition{iCond}.RequiredFilled)
                    if ~isfield(jsonOut,bidsPar.ASLCondition{iCond}.RequiredFilled{iField}) || isempty(jsonOut.(bidsPar.ASLCondition{iCond}.RequiredFilled{iField}))
                        if isempty(strReportFilled)
                            strReportFilled = bidsPar.ASLCondition{iCond}.RequiredFilled{iField};
                        else
                            strReportFilled = [strReportFilled ', ' bidsPar.ASLCondition{iCond}.RequiredFilled{iField}];
                        end
                    end
                end
                
                % Check the required empty fields
                strReportEmpty = '';
                for iField = 1:length(bidsPar.ASLCondition{iCond}.RequiredEmpty)
                    if isfield(jsonOut,bidsPar.ASLCondition{iCond}.RequiredEmpty{iField}) && ~isempty(jsonOut.(bidsPar.ASLCondition{iCond}.RequiredEmpty{iField}))
                        if isempty(strReportEmpty)
                            strReportEmpty = bidsPar.ASLCondition{iCond}.RequiredEmpty{iField};
                        else
                            strReportEmpty = [strReportEmpty ', ' bidsPar.ASLCondition{iCond}.RequiredEmpty{iField}];
                        end
                    end
                end
                
                % Check the recommended filled fields
                strReportRecommended = '';
                for iField = 1:length(bidsPar.ASLCondition{iCond}.RecommendedFilled)
                    if ~isfield(jsonOut,bidsPar.ASLCondition{iCond}.RecommendedFilled{iField}) || isempty(jsonOut.(bidsPar.ASLCondition{iCond}.RecommendedFilled{iField}))
                        if isempty(strReportRecommended)
                            strReportRecommended = bidsPar.ASLCondition{iCond}.RecommendedFilled{iField};
                        else
                            strReportRecommended = [strReportRecommended ', ' bidsPar.ASLCondition{iCond}.RecommendedFilled{iField}];
                        end
                    end
                end
                
                % One of the dependencies was not fulfilled
                if ~isempty(strReportFilled) || ~isempty(strReportEmpty) || ~isempty(strReportRecommended)
                    % Report the conditional field
                    if isempty(bidsPar.ASLCondition{iCond}.value)
                        fprintf('The field %s is empty, please check the dependencies below:\n',bidsPar.ASLCondition{iCond}.field);
                    elseif islogical(bidsPar.ASLCondition{iCond}.value)
                        if bidsPar.ASLCondition{iCond}.value
                            fprintf('The field %s is true, please check the dependencies below:\n',bidsPar.ASLCondition{iCond}.field);
                        else
                            fprintf('The field %s is false, please check the dependencies below:\n',bidsPar.ASLCondition{iCond}.field);
                        end
                    else
                        fprintf('The field %s is %s, please check the dependencies below:\n',bidsPar.ASLCondition{iCond}.field,bidsPar.ASLCondition{iCond}.value);
                    end
                    
                    % Report the incorrect dependencies
                    if ~isempty(strReportFilled)
                        fprintf('The required fields are missing: %s\n',strReportFilled);
                    end
                    if ~isempty(strReportEmpty)
                        fprintf('The following fields should be empty: %s\n',strReportEmpty);
                    end
                    if ~isempty(strReportRecommended)
                        fprintf('The recommended fields are missing: %s\n',strReportRecommended);
                    end
                    fprintf('\n');
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort fields in a predefined order

% Create the structure with the correct field order
fieldOrderStruct = [];
for iField=1:length(bidsPar.listFieldOrder)
    fieldOrderStruct.(bidsPar.listFieldOrder{iField}) = '';
end

% And sort the fields
jsonOut = xASL_adm_OrderFields(jsonOut,fieldOrderStruct);
end
