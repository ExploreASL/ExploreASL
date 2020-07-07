function outStruct = xASL_adm_OrderFields(inStruct,orderStruct)
%xASL_adm_OrderFields Order fields in the structure inStruct to match orderStruct, unmatching fields in inStruct are copied as they are
% at the end, unmatching fields in orderStruct are ignored. This is just a cosmetic change and no values are edited
%
% FORMAT:       outStruct = xASL_adm_OrderFields(inStruct,orderStruct)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Order fields in the structure {{inStruct}} to match
%               {{orderStruct}}, unmatching fields in inStruct are copied as
%               they are at the end, unmatching fields in {{orderStruct}} are 
%               ignored. This is just a cosmetic change and no values are
%               edited.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL


orderArray = fieldnames(orderStruct);
inArray = fieldnames(inStruct);
outStruct ='';
for ii = 1:length(orderArray)
	if isfield(inStruct,orderArray{ii})
		outStruct.(orderArray{ii}) = inStruct.(orderArray{ii});
	end
end

for ii = 1:length(inArray)
	if ~isfield(outStruct,inArray{ii})
		outStruct.(inArray{ii}) = inStruct.(inArray{ii});
	end
end
end
