function json = xASL_io_ReadJson(pathJSON)
%xASL_io_ReadJson Reads the content of a JSON file and converts it to a Matlab structure
%
% FORMAT: json = xASL_io_ReadJson(pathJSON)
%
% INPUT:
%   pathJSON           - path to the JSON file (REQUIRED)
%
% OUTPUT:
%   json               - output Matlab structure with the decoded JSON
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function reads a JSON file located at the given path and decodes it using the Matlab routine
%              JSONDECODE into a Matlab structure.
%
%
% EXAMPLE: json = xASL_io_ReadJson('/tmp/test.json');
% __________________________________
% Copyright (C) 2015-2023 ExploreASL


%% ------------------------------------------------------------------------------------------------
%% 0.   Administration
if nargin < 1 || isempty(pathJSON)
	error('Requires the path to the JSON file.');
end

if ~xASL_exist(pathJSON,'file')
	error(['Cannot find the JSON file: ' pathJSON]);
end

%% 1. Read the JSON file
fileID = fopen(pathJSON, 'r');
try 
	txt = fscanf(fileID, '%c');

	%% 2. Decode the text content
	json = jsondecode(txt);
catch ME
	% In case an error appears during reading or decoding, we close the file and exit
	fclose(fileID);
	error('%s',ME.getReport());
end

% Close the file and returns
fclose(fileID);

end
