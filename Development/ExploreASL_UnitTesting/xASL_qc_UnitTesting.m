function RESULTS = xASL_qc_UnitTesting
%xASL_qc_UnitTesting Script to run all module and submodule tests
%
% FORMAT:       RESULTS = xASL_qc_UnitTesting
% 
% INPUT:        None
%
% OUTPUT:       RESULTS structure
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Script to run all module and submodule tests. Please start
%               this script from the ExploreASL directory.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLES:     RESULTS = xASL_qc_UnitTesting;
% __________________________________
% Copyright 2015-2020 ExploreASL

% Improve command window output
BreakString = [repmat('=',1,100),'\n'];

%% Update GIT
try
    fprintf('GIT: ');
    system('git fetch','-echo');
    system('git pull','-echo');
catch
    warning('It seems that your working directory is not the ExploreASL directory...');
end


%% GET DIRECTORIES
test_parameter_file = fullfile(pwd,'Development','ExploreASL_UnitTesting','xASL_test_parameters.json');
if xASL_exist(test_parameter_file)
    val = jsondecode(fileread(test_parameter_file));
else
    % If the file does not exist yet, it needs to be created
    val = struct;
end

% Check if the current directory is the ASL directory
if xASL_exist(fullfile(pwd,'ExploreASL_Master'))
    val.xASLdir = 'pwd';
else
    val.xASLdir = uigetdir(pwd, 'Select ExploreASL directory...');
end

% Define test directory
val.testDir = uigetdir(pwd, 'Select testing directory...');

% Write the directories to the JSON file
valNew = jsonencode(val);
fid = fopen(test_parameter_file, 'w');
fprintf(fid, '%s', valNew);
fclose(fid);


%% RUN TESTS: INITIALIZATION
fprintf(BreakString);
RESULTS.INIT = runtests('xASL_Initialize_test');
fprintf(BreakString);


%% RUN TESTS: MODULES
fprintf(BreakString);
RESULTS.MODULES.STRUCTURAL = runtests('xASL_module_Structural_test');
fprintf(BreakString);

% RESULTS.MODULES.ASL = runtests('xASL_module_ASL_test');
% RESULTS.MODULES.POPULATION = runtests('xASL_module_Population_test');


%% RUN TESTS: SUBMODULES
% fprintf(BreakString);
% RESULTS.SUBMODULES.xASL_wrp_CreateAnalysisMask = runtests('xASL_wrp_CreateAnalysisMask_test');
% fprintf(BreakString);






end


