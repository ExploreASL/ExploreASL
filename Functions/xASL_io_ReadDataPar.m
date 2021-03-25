function [x] = xASL_io_ReadDataPar(DataParFile)
% xASL_io_ReadDataPar This function reads in a DATA_PAR file and creates the x structure.
%
% FORMAT:   [x] = xASL_io_ReadDataPar(DataParFile)
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
% EXAMPLE:      xASL_io_ReadDataPar('DataParFile.json')
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
% EXAMPLE 3: To define cell arrays, use the official JSON array notation:
%
%   ...
%   "S":
%	{"Atlases": ["TotalGM","DeepWM","Hammers","HOcort_CONN","HOsub_CONN","Mindboggle_OASIS_DKT31_CMA"]}
%   ...
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2021 ExploreASL

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
    % Check if the current field is valid (char, numeric value, structure or cell array)
    if ~(ischar(x.(sFields{n})) || isstruct(x.(sFields{n})) || isnumeric(x.(sFields{n})) || islogical(x.(sFields{n})) || iscell(x.(sFields{n})))
        fprintf('\n%s\n',char(sFields{n}));
        warning('JSON field type could be invalid...');
    end
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

end

%% ConvertNumericalFields (Convert strings which contain numbers to numbers, remove invalid fields)
function [StructOut] = ConvertNumericalFields(StructIn)   

    % Get fields
    listFields = fields(StructIn);
    
    % Iterate over fields
    for iField=1:length(listFields)
        
		% Convert string to num
		TempField = xASL_str2num(StructIn.(listFields{iField}));
            
		% Check if we got an imaginary number
		if max(isnumeric(TempField))
			if imag(TempField)~=0
				TempField = StructIn.(listFields{iField}); % Do not convert complex numbers
			end
		end
            
		% Check if conversion was correct
		if min(isnumeric(TempField)) && size(TempField,1)>size(TempField,2)
			% Enforce a horizontal vector if numerical
			TempField = TempField';
		end

		% Check if there is a numerical vector inside the string
		if min(isnumeric(TempField)) && min(~isnan(TempField))
			StructOut.(listFields{iField}) = TempField;
		else
			% Otherwise keep the string
			StructOut.(listFields{iField}) = StructIn.(listFields{iField});
		end
    end

end
