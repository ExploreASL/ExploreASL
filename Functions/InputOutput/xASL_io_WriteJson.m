function xASL_io_WriteJson(pathJSON, json, bOverwrite)
%xASL_io_WriteJson Writes the content of a json structure to a JSON file
%
% FORMAT: xASL_io_WriteJson(pathJSON, json[, bOverwrite])
%
% INPUT:
%   pathJSON     - path to the JSON file (STRING, REQUIRED)
%   json         - Matlab structure to be written as JSON (STRUCT, REQUIRED)
%   bOverwrite   - boolean that is true if an existing file should be overwritten (BOOLEAN, OPTIONAL, DEFAULT = TRUE)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function writes a Matlab structure into a JSON file located at the given path encoding it using the Matlab function
%              JSONENCODE.
%
%
% EXAMPLE: xASL_io_WriteJson('/tmp/test.json', json);
% __________________________________
% Copyright (C) 2015-2023 ExploreASL


%% ------------------------------------------------------------------------------------------------
%% 0.   Administration
if nargin < 1 || isempty(pathJSON)
	error('The path to the JSON file is missing.');
end

if nargin < 2 || isempty(json)
	error('The json structure to be saved is missing.')
end

if nargin < 3 || isempty(bOverwrite)
	bOverwrite = true;
end

if xASL_exist(pathJSON,'file') && ~bOverwrite
	error(['File exists, please set bOverwrite to true: ' pathJSON]);
end

%% 1. Encode the content to text.
% This is executed first in case there are issues with the file
% Our version of PrettyPrint is preferred to avoid printing each vector to several lines
%if verLessThan('matlab', '9.10')
	txt = jsonencode(json);
	txt = xASL_sub_PrettyPrint(txt);
%else
%	txt = jsonencode(json, 'PrettyPrint', true);
%end

%% 2. Open the JSON file for writing and save it
fileID = fopen(pathJSON,'w+');

try
	fprintf(fileID, '%s', txt);
catch ME
	% In case an error appears during writing, we close the file and exit
	fclose(fileID);
	error('%s',ME.getReport());
end

% Close the file and return
fclose(fileID);

end


%% Add newlines to make json output readable
function txt = xASL_sub_PrettyPrint(txt)

	txt = strrep(txt, '",', ['",'  newline]);
	txt = strrep(txt, '],', ['],'  newline]);
	txt = strrep(txt, '},', ['},'  newline]);

end

