function [x] = xASL_imp_BIDS2Legacy(x)
%xASL_imp_BIDS2Legacy BIDS2LEGACY conversion script which calls xASL_bids_BIDS2Legacy.
%
% FORMAT: [x] = xASL_imp_BIDS2Legacy(x);
%
% INPUT:
%   x             - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x             - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    BIDS to Legacy conversion script which calls xASL_bids_BIDS2Legacy.
%
% 1. Start with checking dataset_description.json & rawdata
% - 1. The input is dataset_description.json in the rawdata folder
% - 2. The input is dataPar.json or sourceStructure.json - have to look for a rawdata folder
% 2. Run the legacy conversion: Check if a dataPar is provided, otherwise use the defaults
% 3. Overwrite DataParPath
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% 1. Start with checking dataset_description.json & rawdata
    % Either the rawdata/dataset_description.json was provided
	% Or the dataPar.json was provided and we have to search was dataset_description.json

	[Fpath, Ffile, Fext] = fileparts(x.DataParPath);
	if strcmp([Ffile Fext], 'dataset_description.json')
		%% 1.1 The input is dataset_description.json in the rawdata folder
		[localStudyRoot, Ffile] = fileparts(Fpath);
		if ~strcmp(Ffile, 'rawdata')
			error('Invalid folder in which dataset_description.json was found, should be /rawdata');
		end
	else
		%% 1.2 The input is dataPar.json or sourceStructure.json - have to look for a rawdata folder
		pathRawData = fullfile(x.dir.StudyRoot,'rawdata');
		% Check the if the correct BIDS structure is in place
		if ~exist(pathRawData,'dir')
			error('Path rawdata does not exist');
		elseif ~exist(fullfile(pathRawData,'dataset_description.json'),'file')
			error('File dataset_description.json is not found');
		else
			localStudyRoot = x.dir.StudyRoot;
		end
	end

	%% 2. Run the legacy conversion: Check if a dataPar is provided, otherwise use the defaults
	fListDataPar = xASL_adm_GetFileList(localStudyRoot,'(?i)(^dataPar.*\.json$)', 'FPList', [], 0);
	if length(fListDataPar) < 1
		fprintf('There is no dataPar.json file in the study root directory. Default settings will be used...\n');
		% Fill the dataPars with default parameters
		dataPar = xASL_bids_BIDS2Legacy(localStudyRoot, 1, []);
	else
		% Fill the dataPars with the provided parameters
		dataPar = spm_jsonread(fListDataPar{1});
		dataPar = xASL_bids_BIDS2Legacy(localStudyRoot, 1, dataPar);
	end

	%% 3. Overwrite DataParPath
	x.DataParPath = dataPar.x.DataParPath;
    
    
end


