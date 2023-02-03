function xASL_io_WriteJson(pathJSON, json, bOverwrite)
%xASL_io_WriteJson Writes the content of a json structure to a JSON file
%
% FORMAT: xASL_io_WriteJson(pathJSON, json[, bOverwrite])
%
% INPUT:
%   pathJSON     - path to the JSON file (REQUIRED)
%   json         - Matlab structure to be written as JSON (REQUIRED)
%   bOverwrite   - boolean that is true if an existing file should be overwritten (OPTIONAL, DEFAULT = TRUE)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function writes a Matlab structure into a JSON file located at the given path encoding it using the Matlab routine
%              JSONENCODE.
%
%
% EXAMPLE: xASL_io_WriteJson('/tmp/test.json', json);
% __________________________________
% Copyright (C) 2015-2023 ExploreASL


%% ------------------------------------------------------------------------------------------------
%% 0.   Administration
if nargin < 1 || isempty(pathJSON)
	error('Requires the path to the JSON file.');
end

if nargin < 2 || isempty(json)
	error('Requires the json structure to be saved.')
end

if nargin < 3 || isempty(bOverwrite)
	bOverwrite = true;
end

if xASL_exist(pathJSON,'file') && ~bOverwrite
	error(['File exists, please set bOverwrite to true: ' pathJSON]);
end

%% 1. Open the JSON file for writing
fileID = fopen(pathJSON,'w+');


%% 2. Encode the content to text and write it to the file
txt = jsonencode(json);
fprintf(fileID, '%s', txt);

% Close the file and return
fclose(fileID);

end
