function [json] = xASL_bids_CreateDatasetDescriptionTemplate(draft, versionExploreASL)
%xASL_bids_CreateDatasetDescriptionTemplate This script creates a JSON structure which can be saved
% using spm_jsonwrite to get a dataset_description.json template.
%
% FORMAT: [json] = xASL_bids_CreateDatasetDescriptionTemplate(draft)
% 
% INPUT:
%   draft             - Structure which defines the dataset_description fields (STRUCT, REQUIRED)
%   versionExploreASL - ExploreASL version (STRING, REQUIRED)
%
% OUTPUT:
%   json        - JSON structure which can be saved using spm_jsonwrite (dataset_description.json template)
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This script creates a JSON structure which can be saved
%               using spm_jsonwrite to get a dataset_description.json template.
%               Missing fields that are required are added. BIDSVersion checked against the current configured version.
%               Remaining fields will be validated. Other fields not belonging to dataset_description.json are ignored.
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

    % Create the output dataset_description.json
    json = struct;
    
	% Add missing required fields
    for iCell = 1:numel(bidsPar.datasetDescription.Required)
		% Checks if the provided version field matches with the currently configured version
		if strcmp(bidsPar.datasetDescription.Required{1,iCell},'BIDSVersion')
			if isfield(draft,'BIDSVersion')
				if ~strcmp(draft.BIDSVersion,bidsPar.BIDSVersion)
					warning('Difference between the provided and current BIDSVersion');
				end
				json.BIDSVersion = draft.BIDSVersion;
			else
				json.BIDSVersion = bidsPar.BIDSVersion;
			end
		else
			% Other required fields - just copy or label as Undefined
			if isfield(draft,bidsPar.datasetDescription.Required{1,iCell})
				json.(bidsPar.datasetDescription.Required{1,iCell}) = draft.(bidsPar.datasetDescription.Required{1,iCell});
			else
				fprintf('Warning: Add undefined required field %s\n',bidsPar.datasetDescription.Required{1,iCell});
				json.(bidsPar.datasetDescription.Required{1,iCell}) = 'Undefined';
			end
		end
	end
	
	listMissingFiles = [];
	
	% Recommended fields
	for iCell = 1:numel(bidsPar.datasetDescription.Recommended)
		% Only copy if defined
		if isfield(draft,bidsPar.datasetDescription.Recommended{1,iCell})
			json.(bidsPar.datasetDescription.Recommended{1,iCell}) = draft.(bidsPar.datasetDescription.Recommended{1,iCell});
		else
			if length(listMissingFiles)>1
				listMissingFiles = [listMissingFiles ', '];
			end
			listMissingFiles = [listMissingFiles bidsPar.datasetDescription.Recommended{1,iCell}];
		end
	end
	
	if length(listMissingFiles)>1
		% Report the missing fields
        fprintf('================================== dataset_description.json ==================================\n');
		fprintf('Missing recommended fields:           %s\n',listMissingFiles);
	end
		
    % Optional fields
    for iCell = 1:numel(bidsPar.datasetDescription.Optional)
        % Only copy if defined
        if isfield(draft,bidsPar.datasetDescription.Optional{1,iCell})
            json.(bidsPar.datasetDescription.Optional{1,iCell}) = draft.(bidsPar.datasetDescription.Optional{1,iCell});
        end
        % Add version which was used for the import to the Acknowledge field (Imported with xASL v1.x.x)
        if strcmp(bidsPar.datasetDescription.Optional{1,iCell},'Acknowledgements')
            json.Acknowledgements = ['Imported with xASL ' versionExploreASL];
        end
    end
    
    % Add linebreak after printing
    if length(listMissingFiles)>1
        fprintf('\n');
    end
    
end
