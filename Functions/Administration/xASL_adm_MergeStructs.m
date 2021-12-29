function mergedStruct = xASL_adm_MergeStructs(mainStruct, secondaryStruct)
% xASL_adm_MergeStructs Adds all fields from the secondaryStruct to the mainStruct while avoid overwriting
%
% FORMAT: mergedStruct = xASL_adm_MergeStructs(mainStruct, secondaryStruct)
%
% INPUT:
%   mainStruct      - the structure with master information
%   secondaryStruct - the structure that is merged in
%
% OUTPUT:
%   mergedStructure - the final merged structure
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: It merges two structures. It takes everything from the mainStruct and keep it as it is. It adds all fields from the secondaryStructure
%              to the main structure while checking for duplicates. It is not overwriting anything, all duplicit content is taken from mainStruct.
%              It works iteratively by correctly merging also the substructs.
%
% EXAMPLE:     mergedStruct = xASL_adm_MergeStructs(mainStruct, secondaryStruct)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:  n/a
%
% Copyright 2015-2021 ExploreASL

% Take the mainStruct as the basis
mergedStruct = mainStruct;

fieldsList = fields(secondaryStruct);
for iField=1:length(fieldsList)
	if ~isfield(mainStruct,(fieldsList{iField}))
		% Field was missing in the main, so it can be taken from the secondary
		mergedStruct.(fieldsList{iField}) = secondaryStruct.(fieldsList{iField});
	elseif isstruct(mainStruct.(fieldsList{iField})) && isstruct(secondaryStruct.(fieldsList{iField}))
		% The field is in both and it is a structure - then re-run this function iterative on the substruct
		mergedStruct.(fieldsList{iField}) = xASL_adm_MergeStructs(mainStruct.(fieldsList{iField}),secondaryStruct.(fieldsList{iField}));
	end
	% In case both structs have the field and the field is not a substract then keep the content from mainStruct and do not do anything
end
