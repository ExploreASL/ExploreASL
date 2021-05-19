function [TestResults] = xASL_test_IntegrationTesting
%xASL_test_IntegrationTesting Main integration testing script
%
% INPUT:        n/a
%
% OUTPUT:       TestResults  - Test results struct
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Integration testing script. Run the full pipeline on the flavor datasets using ExploreASL_Master.
%
% EXAMPLE:      [TestResults] = xASL_test_IntegrationTesting;
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


    %% Initialize Paths

    % Add testing directory
    TestingPath = mfilename('fullpath');    % Path of the current script
    filepath = fileparts(TestingPath);      % Testing directory
    xASL_dir = fileparts(filepath);         % All subfolders of xASL directory
    warning('off','all')                    % Turn of warnings (only add path related)
    addpath(genpath(xASL_dir))              % Add all scripts to path (including testing directory)
    warning('on','all')

    %% Clean Up
    clc
    clear

    %% Get paths
    if exist(fullfile(TestingPath,'testConfig.json'),'file')
    	testConfig = spm_jsonread(fullfile(TestingPath,'testConfig.json'));
    else
    	fprintf('Please add a testConfig.json to the Testing directory of ExploreASL...\n');
    	return
    end

    % Check testConfig fields
    if ~isfield(testConfig, 'pathExploreASL')
    	error('Add the pathExploreASL field to your testConfig.json...');
    end
    if ~isfield(testConfig, 'pathTest')
    	error('Add the pathTest field to your testConfig.json...');
    end
    if ~isfield(testConfig, 'cmdCloneFlavors')
    	error('Add the cmdCloneFlavors field to your testConfig.json...');
    end

    %% Set-up
	cd(testConfig.pathExploreASL);

	% Initialize ExploreASL
	[x] = ExploreASL_Initialize;
	cd(pathTest);
	
	% Clone flavors if they do not exist already
	if ~exist(fullfile(testConfig.pathTest,'FlavorDatabase'), 'dir')
	    xASL_adm_CreateDir(testConfig.pathTest);
	    system(testConfig.cmdCloneFlavors);
	end

	% Get the up-to-date flavor list
	% ...

	%% Iterate over flavors
	% ...

	% Run individual flavor
	% ...

    

end



