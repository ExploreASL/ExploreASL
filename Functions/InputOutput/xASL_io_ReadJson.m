function json = xASL_io_ReadJson(pathJSON)
%xASL_io_ReadJson Reads the content of a JSON file and converts it to a Matlab structure
%
% FORMAT: json = xASL_io_ReadJson(pathJSON)
%
% INPUT:
%   pathJSON           - path to the JSON file (STRING, REQUIRED)
%
% OUTPUT:
%   json               - output Matlab structure with the decoded JSON
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function reads a JSON file located at the given path and decodes it using the Matlab built-in function JSONDECODE. The output is a Matlab
% structure with the decoded JSON content.
%
% EXAMPLE: json = xASL_io_ReadJson('/tmp/test.json');
% __________________________________
% Copyright (C) 2015-2023 ExploreASL


%% ------------------------------------------------------------------------------------------------
%% 0.   Administration
if nargin < 1 || isempty(pathJSON)
	error('The path to the JSON file is missing.'); % If the path is missing, we exit with an error
end

if ~xASL_exist(pathJSON,'file')
	error(['Cannot find the JSON file: ' pathJSON]);
end

%% ------------------------------------------------------------------------------------------------
%% 1. Load the JSON file
fileID = fopen(pathJSON, 'r'); % Open the file for reading

try 
	txt = fscanf(fileID, '%c'); % Read the content of the file
catch ME
	% In case an error appears during reading or decoding, we close the file and exit
	fclose(fileID); % Close the file
	error('%s',ME.getReport()); % Display the error message
end

% Close the file
fclose(fileID); % Close the file

%% ------------------------------------------------------------------------------------------------
%% 2. Decode the text content

try
	json = jsondecode(txt); % Decode the JSON content
catch ME
	% In case an error appears during reading or decoding, we exit
	fprintf('Something went wrong during the decoding of the JSON file: %s\n', pathJSON');
	fprintf('The JSON file might be corrupted. See the error message below or check the JSON file content at https://jsonlint.com \n');
	error('%s',ME.getReport()); % Display the error message
end

% In case we have a multi-context import and there is equivalent information in all contexts, it will decode each context as a struct instead of a cell
% This needs to be fixed "semi-manually"
if isfield(json, 'StudyPars') && isstruct(json.StudyPars)
	% We need to convert StudyPars field to a cell
	jsonOld = json; % Save the original JSON
	json = rmfield(json, 'StudyPars'); % Keep all fields but StudyPars in the new version
	for iField = 1:numel(jsonOld.StudyPars)
		json.StudyPars{iField,1} = jsonOld.StudyPars(iField);
	end
end

end