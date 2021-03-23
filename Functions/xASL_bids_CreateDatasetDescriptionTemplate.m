function [json] = xASL_bids_CreateDatasetDescriptionTemplate(varargin)
%xASL_bids_CreateDatasetDescriptionTemplate This script creates a JSON structure which can be saved
% using spm_jsonwrite to get a dataset_description.json template.
%
% FORMAT: [json] = xASL_bids_CreateDatasetDescriptionTemplate(varargin)
% 
% INPUT:
%   Name                - Patient name (CHAR ARRAY, OPTIONAL, DEFAULT = 'Undefined')
%   BIDSVersion         - BIDS Version (e.g. 1.5.0) (CHAR ARRAY, OPTIONAL, DEFAULT = use xASL_bids_Config to get version)
%   HEDVersion          - HED Version (CHAR ARRAY, OPTIONAL, DEFAULT = 'Undefined')
%   DatasetType         - Dataset Type (CHAR ARRAY, OPTIONAL, DEFAULT = 'Undefined')
%   License             - License (CHAR ARRAY, OPTIONAL, DEFAULT = 'Undefined')
%   Authors             - Authors (CHAR ARRAY, OPTIONAL, DEFAULT = 'Undefined')
%   Acknowledgements    - Acknowlegements (CHAR ARRAY, OPTIONAL, DEFAULT = 'Undefined')
%
% OUTPUT:
%   json                - JSON structure which can be saved using spm_jsonwrite (dataset_description.json template)
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This script creates a JSON structure which can be saved
%               using spm_jsonwrite to get a dataset_description.json template.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      [json] = xASL_bids_CreateDatasetDescriptionTemplate('Test Patient','1.5.0','1.2.3','Test','Test License','Author Names','Thanks');
%               [json] = xASL_bids_CreateDatasetDescriptionTemplate();
%               
% __________________________________
% Copyright 2015-2021 ExploreASL

    % Get default BIDS configuration
    bidsPar = xASL_bids_Config();

    % Create dummy dataset_description.json
    json = struct;
    
    % Required fields
    for iCell = 1:numel(bidsPar.datasetDescription.Required)
         json.(bidsPar.datasetDescription.Required{1,iCell}) = 'Undefined';
    end
    % Recommended fields
    for iCell = 1:numel(bidsPar.datasetDescription.Recommended)
         json.(bidsPar.datasetDescription.Recommended{1,iCell}) = 'Undefined';
    end
    % Optional fields
    for iCell = 1:numel(bidsPar.datasetDescription.Required)
         json.(bidsPar.datasetDescription.Optional{1,iCell}) = 'Undefined';
    end
    
    % Parse the input parameters
    p = inputParsing(varargin{:});
    
    % Optionally change the parameter values
    json.Name = p.Results.Name;
    json.BIDSVersion = p.Results.BIDSVersion;
    json.HEDVersion = p.Results.HEDVersion;
    json.DatasetType = p.Results.DatasetType;
    json.License = p.Results.License;
    json.Authors = p.Results.Authors;
    json.Acknowledgements = p.Results.Acknowledgements;
    
end

%% Parse input
function p = inputParsing(varargin)

    % Get default BIDS configuration
    bidsPar = xASL_bids_Config();

    % Initialize input parser
    p = inputParser;
    
    % Define valid input variables
    validName = @(variable) ischar(variable);
    validBIDSVersion = @(variable) ischar(variable);
    validHEDVersion = @(variable) ischar(variable);
    validDatasetType = @(variable) ischar(variable);
    validLicense = @(variable) ischar(variable);
    validAuthors = @(variable) ischar(variable);
    validAcknowledgements = @(variable) ischar(variable);
    
    % Define defaults
    defaultName = 'Undefined';
    defaultBIDSVersion = bidsPar.BIDSVersion;
    defaultHEDVersion = 'Undefined';
    defaultDatasetType = 'Undefined';
    defaultLicense = 'Undefined';
    defaultAuthors = 'Undefined';
    defaultAcknowledgements = 'Undefined';
    
    % Add definitions to the input parser
    addOptional(p, 'Name', defaultName, validName);
    addOptional(p, 'BIDSVersion', defaultBIDSVersion, validBIDSVersion);
    addOptional(p, 'HEDVersion', defaultHEDVersion, validHEDVersion);
    addOptional(p, 'DatasetType', defaultDatasetType, validDatasetType);
    addOptional(p, 'License', defaultLicense, validLicense);
    addOptional(p, 'Authors', defaultAuthors, validAuthors);
    addOptional(p, 'Acknowledgements', defaultAcknowledgements, validAcknowledgements);
    
    % Parse input
    parse(p,varargin{:});

end



