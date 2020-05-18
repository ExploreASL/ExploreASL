function [PathJSON] = xASL_init_ConvertM2JSON(PathM, bOverwrite)
%xASL_init_ConvertM2JSON Convert the old DataPar m-file to JSON format
%
% FORMAT: [PathJSON] = xASL_init_ConvertM2JSON(PathM)
%
% INPUT:
%   PathM       - path to legacy data parameterfile of ExploreASL (REQUIRED)
%   bOverwrite  - boolean specifying if existing PathJSON will be
%                 overwritten or not (OPTIONAL, DEFAULT=true)
%
% OUTPUT:       - path to JSON file containing parameters for running
%                 ExploreASL
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function converts and replaces the legacy data parameter m-format
%              by a JSON file. A DataPar.m was the settings/parameter file, specific to a dataset to be
%              processed by ExploreASL, now replaced to JSON by BIDS.
%              Note that the deployed/compiled version of ExploreASL
%              requires the JSON file, this function should not be compiled
%              along. This function performs the following steps:
%
%              1) Run the m-file to load parameters
%              2) Escape characters that are illegal in JSON
%              3) Write the JSON
%
% EXAMPLE: PathJSON = xASL_init_ConvertM2JSON('/MyStudy/DataParameterFile.m');
% __________________________________
% Copyright (C) 2015-20120 ExploreASL


%% -----------------------------------------------
%% Admin
if nargin<2 || isempty(bOverwrite)
    bOverwrite = true;
end

[Fpath, Ffile] = fileparts(PathM);
CurrentFolder = pwd;
if ~isempty(Fpath)
    cd(Fpath);
end
PathJSON = fullfile(Fpath, [Ffile '.json']);

%% -----------------------------------------------
%% 1) Run the m-file to load parameters
FunctionHandle = str2func(Ffile);
x = FunctionHandle();

cd(CurrentFolder);

%% -----------------------------------------------
%% 2) Escape characters that are illegal in JSON
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
        
        % Escape illegal characters
		if ~iscell(tString)
			Strs = find(tString=='\');
			for iStr=1:length(Strs)
				tString = [tString(1:Strs(iStr)-1) '\\' tString(Strs(iStr)+1:end)];
				Strs = Strs+1;
			end
		else
			for ii=1:length(tString)
				Strs = find(tString{ii}=='\');
				for iStr=1:length(Strs)
					tString{ii} = [tString{ii}(1:Strs(iStr)-1) '\\' tString{ii}(Strs(iStr)+1:end)];
					Strs = Strs+1;
				end
			end
		end
        x.(FieldsX{iField}) = tString;
    end
end

%% -----------------------------------------------
%% 3) Write the JSON
if bOverwrite
    xASL_delete(PathJSON); % remove any previous version
elseif ~bOverwrite && exist(PathJSON, 'file')
    warning([PathJSON ' already existed, skipping']);
    return;
end

spm_jsonwrite(PathJSON, x);

end

