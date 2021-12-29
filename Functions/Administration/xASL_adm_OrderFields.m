function outStruct = xASL_adm_OrderFields(inStruct,orderStruct)
%xASL_adm_OrderFields Order fields in the structure inStruct to match orderStruct, unmatching fields in inStruct are copied as they are
% at the end, unmatching fields in orderStruct are ignored. This is just a cosmetic change and no values are edited
%
% FORMAT:       outStruct = xASL_adm_OrderFields(inStruct,orderStruct)
% 
% INPUT:        inStruct      - Struct that should be re-ordered (STRUCT, REQUIRED)
%               orderStruct   - Struct defining the order (STRUCT, REQUIRED)
%
% OUTPUT:       outStruct     - output struct (STRUCT)
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Order fields in the structure {{inStruct}} to match {{orderStruct}},
%               unmatching fields in inStruct are copied as
%               they are at the end, unmatching fields in {{orderStruct}} are 
%               ignored. This is just a cosmetic change and no values are
%               edited.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:    
%
% struct1.A = 'ABC';
% struct1.B = 'DEF';
% struct1.C = 'GHI';
% struct1.D = 'JKL';
% struct2.A = 1;
% struct2.B = 2;
% struct2.C = 3;
% outStruct = xASL_adm_OrderFields(struct1,struct2);
% clear struct2
% 
% struct2.C = 3;
% struct2.B = 2;
% struct2.A = 1;
% outStruct = xASL_adm_OrderFields(struct1,struct2);
% __________________________________
% Copyright 2015-2021 ExploreASL


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
