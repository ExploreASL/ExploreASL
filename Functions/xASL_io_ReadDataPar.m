function [x] = xASL_io_ReadDataPar(DataParFile)
% xASL_io_ReadDataPar This function reads DATA_PAR .json or .m file and creates the x structure.
%
% FORMAT:   [x] = xASL_io_ReadDataPar(DataParFile)
%
% INPUT:
%   DataParFile     - Filename of the DATA_PAR file with either .m or .json extension
%
% OUTPUT:
%   x               - x structure
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function reads in a DATA_PAR .json or .m file and creates the x
%               structure. The name of the DATA_PAR file is given as a string or
%               character array. The output is the x structure. It only loads the data, removes the x-prefixes, 
%               but keeps all the field names and units. It doesn't do any conversions to or from BIDS. The 
%               only added value to normal json-read is that it detects invalid entries (numbers in strings, and 
%               weird arrays), converts them correctly and reports this incorrect entries so that they can be manually 
%               fixed. Also, if an .m file is provided, it converts and saves it to a JSON file (doesn't overwrite)
%               and reports that you should stop using .m files.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:      xASL_io_ReadDataPar('DataParFile.json')
%
% EXAMPLE 2:    xASL_io_ReadDataPar('DataParFile.m')
%
% EXAMPLE 2:    JSON FILE
%
% {
% 	"x": [{
% 			"name": 				"ExampleDataSet",
% 			"subject_regexp": 		"^Sub-\\d{3}$",
% 			"M0": 					"separate_scan",
% 			"BackgroundSuppressionNumberPulses": 2,
% 			"readout_dim": 			"2D",
% 			"QUALITY": 				false,
% 			"Vendor": 				"Philips",
% 		}]
% }
%
% Don't forget to escape the backslashes!
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
if nargin < 1 || isempty(DataParFile) || ~exist(DataParFile, 'file')
    error('DataParFile does not exist...');
end

% Input has to be a character array
DataParFile = char(DataParFile);
[Fpath, Ffile, Fext] = fileparts(DataParFile);

%% Read an m-file
if strcmp(Fext, '.m')
	% Read the m-file by executing it
	CurrentFolder = pwd;
	if ~isempty(Fpath)
		cd(Fpath);
	end

	FunctionHandle = str2func(Ffile);
	x = FunctionHandle();

	cd(CurrentFolder);
	
	PathJSON = fullfile(Fpath, [Ffile '.json']);
	
	% Escape characters that are illegal in JSON
	% Note that most characters (like backslashes) are correctly escaped when saving the JSON in spm_jsonwrite
	FieldsX = fields(x);
	for iField=1:length(FieldsX)
		tString = x.(FieldsX{iField});
		
		if ~isstruct(tString)
			% fix Boolean conversion issue JSON
			if islogical(tString) && tString
				tString = 1;
			elseif islogical(tString) && ~tString
				tString = 0;
			end
			
			x.(FieldsX{iField}) = tString;
		end
	end
	% Write the JSON if it doesn't exist
	if exist(PathJSON, 'file')
		fprintf('Warning: m-files are deprecated, but cannot convert to json-file because a file with json extension exists in the same directory');
	else
		fprintf('Warning: m-files are deprecated, converted and saved to json-file. Please delete the original m-file');
		spm_jsonwrite(PathJSON, x);
	end

elseif strcmp(Fext, '.json')
	%% Read a JSON file	
	jsonData = spm_jsonread(DataParFile);

	% Remove x field if it exists
	if isfield(jsonData,'x')
		x = jsonData.x;
	else
		x = jsonData;
	end
    
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

	% Convert strings containing numbers to number
	x = xASL_io_ReadDataPar_FixFields(x);
else
	error('Unknown file extension');
end

end

%% Goes through the JSON fields and make sure that incorrect entries are corrected
function [StructOut] = xASL_io_ReadDataPar_FixFields(StructIn)

% First fixes the numerical fields in the root structure and x.Q substructure
StructOut = xASL_io_ReadDataPar_ConvertNumericalFields(StructIn);
if isfield(StructIn,'Q')
	StructOut.Q = xASL_io_ReadDataPar_ConvertNumericalFields(StructIn.Q);
end

% Looks for fields with a specific cell array structure entered as a string, which is incorrect by new JSON standards, and fixes this
StructOut = xASL_io_ReadDataPar_ConvertCellArrays(StructOut);

end

%% ConvertNumericalFields (Convert strings which contain numbers to numbers, remove invalid fields)
function [StructOut] = xASL_io_ReadDataPar_ConvertNumericalFields(StructIn)   

% Get fields
listFields = fields(StructIn);

% Iterate over fields
for iField=1:length(listFields)
	
	% Checks if conversion should be attempted. This check is done also inside the xASL_str2num, however,
	% it is good to do it also beforehand, because you avoid unnecessary processing of numeric strings which
	% can result in false-warnings about numeric being given inside a string. 
	% i.e. this check avoids extra processing and warnings eventhough functionally not necessary
	if ~isempty(StructIn.(listFields{iField})) && ~isnumeric(StructIn.(listFields{iField}))
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
			if ~isequal(TempField,StructIn.(listFields{iField}))
				fprintf('Warning: A numerical field provided as a string (%s). Converted, but should be avoided.\n', StructIn.(listFields{iField}));
			end
		else
			% Otherwise keep the string
			StructOut.(listFields{iField}) = StructIn.(listFields{iField});
		end
	else
		StructOut.(listFields{iField}) = StructIn.(listFields{iField});
	end
end

end

%%xASL_io_ReadDataPar_ConvertCellArrays finds all the cell arrays entered as strings and fixes them
function StructOut = xASL_io_ReadDataPar_ConvertCellArrays(StructIn)

% Get fields
listFields = fields(StructIn);

% Iterate over fields
for iField=1:length(listFields)
	
	% Finds the incorrect fields
	if isstruct(StructIn.(listFields{iField}))
		% Calls recursively the function for structs
		StructOut.(listFields{iField}) = xASL_io_ReadDataPar_ConvertCellArrays(StructIn.(listFields{iField}));
		
	elseif ischar(StructIn.(listFields{iField})) && ~isempty(StructIn.(listFields{iField})) && ...
			StructIn.(listFields{iField})(1) == '{' && StructIn.(listFields{iField})(end) == '}'
		% For chars entered in {} brackets, brakes the entries in '' to individuals arrays with a cell array
		iFirst = 2;
		iApostrophe = strfind(StructIn.(listFields{iField}),39);
		StructOut.(listFields{iField}) = cell(1,0);
		iCurrent = 1;
		while(length(iApostrophe) >= iCurrent)
			% Sets the string start to after the next apostrophe
			iFirst = iApostrophe(iCurrent)+1;
			
			% If a ending apostrophe is present to the current one, then write up the string to the cell and move on
			if length(iApostrophe) > iCurrent
				StructOut.(listFields{iField}){1,length(StructOut.(listFields{iField}))+1} = StructIn.(listFields{iField})(iFirst:(iApostrophe(iCurrent+1)-1));
				% Else we just end a let the remains be written
				iFirst = iApostrophe(iCurrent+1)+1;
			end
			iCurrent = iCurrent + 2;
		end
		% If some opened strings remain, write them to a last array field
		if iFirst<=length(StructIn.(listFields{iField}))-1
			StructOut.(listFields{iField}){1,length(StructOut.(listFields{iField}))+1} = StructIn.(listFields{iField})(iFirst:(length(StructIn.(listFields{iField}))-1));
		end
		fprintf('Warning: An incorrect format of a cell array provided (%s). Converted, but should be avoided.\n', StructIn.(listFields{iField}));
	else
		% Otherwise keep the string
		StructOut.(listFields{iField}) = StructIn.(listFields{iField});
	end
end
end