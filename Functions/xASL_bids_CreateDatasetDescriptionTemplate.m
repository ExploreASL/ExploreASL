function [json] = xASL_bids_CreateDatasetDescriptionTemplate(draft)
%xASL_bids_CreateDatasetDescriptionTemplate This script creates a JSON structure which can be saved
% using spm_jsonwrite to get a dataset_description.json template.
%
% FORMAT: [json] = xASL_bids_CreateDatasetDescriptionTemplate(draft)
% 
% INPUT:
%   draft       - Structure which defines the dataset_description fields (STRUCT, REQUIRED)
%
% OUTPUT:
%   json        - JSON structure which can be saved using spm_jsonwrite (dataset_description.json template)
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This script creates a JSON structure which can be saved
%               using spm_jsonwrite to get a dataset_description.json template.
%               Missing fields that are required are added.
%               Remaining fields will be validated.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      draft.Name = 'DRO_Digital_Reference_Object';
%               draft.Test = 'Test';
%               draft.License = 'Test_License';
%               [json] = xASL_bids_CreateDatasetDescriptionTemplate(draft);
%               
% __________________________________
% Copyright 2015-2021 ExploreASL

    % Get default BIDS configuration
    bidsPar = xASL_bids_Config();

    % Create dummy dataset_description.json
    validFields = struct;
    
    % Required fields
    for iCell = 1:numel(bidsPar.datasetDescription.Required)
         validFields.(bidsPar.datasetDescription.Required{1,iCell}) = 'Undefined';
    end
    % Recommended fields
    for iCell = 1:numel(bidsPar.datasetDescription.Recommended)
         validFields.(bidsPar.datasetDescription.Recommended{1,iCell}) = 'Undefined';
    end
    % Optional fields
    for iCell = 1:numel(bidsPar.datasetDescription.Required)
         validFields.(bidsPar.datasetDescription.Optional{1,iCell}) = 'Undefined';
    end
    
    % Check the draft fieldnames
    fieldsDraft = fieldnames(draft);
    
    % Validate the fields
    for iField = 1:size(fieldsDraft,1)
        % Check if fieldname is valid, otherwise remove it
        if ~isfield(validFields,fieldsDraft{iField})
            fprintf('Remove invalid field %s...\n',fieldsDraft{iField});
            draft = rmfield(draft,fieldsDraft{iField});
        end
    end
    
    % Add missing required fields
    for iCell = 1:numel(bidsPar.datasetDescription.Required)
        % Check required field
        if ~isfield(draft,bidsPar.datasetDescription.Required{1,iCell})
            fprintf('Add required field %s...\n',bidsPar.datasetDescription.Required{1,iCell});
            if strcmp(bidsPar.datasetDescription.Required{1,iCell},'BIDSVersion')
                % Add correct BIDS version
                draft.(bidsPar.datasetDescription.Required{1,iCell}) = bidsPar.BIDSVersion;
            else
                % Add undefined field
                draft.(bidsPar.datasetDescription.Required{1,iCell}) = 'Undefined';
            end
        end
    end
    
    % Export valid JSON struct
    json = draft;
    
    
end




