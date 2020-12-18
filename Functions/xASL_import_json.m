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
% 			"BackGrSupprPulses":	"2",
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

% Input has to be a character array (we convert it everytime now, just to be sure)
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
    if ischar(x.(sFields{n})) % Make sure the field is a character array
        if contains(x.(sFields{n}),'{') && contains(x.(sFields{n}),'}') % Automatically detect cell arrays (like the "Atlases" e.g.)
            tempField = x.(sFields{n});
            tempField = strsplit(tempField,',');
            tempField = strrep(strrep(tempField,'{',''),'}','');
            x.(sFields{n}) = strrep(tempField,'''','');
        end
    end
    % Also look for cell arrays in sub structures (like the atlases in the S struct e.g.)
    if isstruct(x.(sFields{n}))
        sFieldsStruct = fieldnames(x.(sFields{n}));
        for m=1:size(sFieldsStruct,1)
            if ischar(x.(sFields{n}).(sFieldsStruct{m})) % Make sure the field is a character array
                if contains(x.(sFields{n}).(sFieldsStruct{m}),'{') && contains(x.(sFields{n}).(sFieldsStruct{m}),'}') % Automatically detect cell arrays (like the "Atlases" e.g.)
                    tempField = x.(sFields{n}).(sFieldsStruct{m});
                    tempField = strsplit(tempField,',');
                    tempField = strrep(strrep(tempField,'{',''),'}','');
                    x.(sFields{n}).(sFieldsStruct{m}) = strrep(tempField,'''','');
                end
            end
        end
    end
end





%% Convert strings containing numbers to number
x = ConvertNumericalFields(x);
if isfield(x,'Q')
	x.Q = ConvertNumericalFields(x.Q);
end

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