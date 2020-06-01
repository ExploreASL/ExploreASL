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
% DESCRIPTION: Script to run all module and submodule tests
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLES:     RESULTS = xASL_qc_UnitTesting
% __________________________________
% Copyright 2015-2020 ExploreASL

%% Update GIT
try
    fprintf('GIT: ');
    system('git fetch','-echo');
    system('git pull','-echo');
catch
    warning('It seems that your working directory is not the ExploreASL directory...');
end


%% RUN TESTS: MODULES
fprintf([repmat('=',1,100),'\n']);
RESULTS.MODULES.STRUCTURAL = runtests('xASL_module_Structural_test');
fprintf([repmat('=',1,100),'\n']);

% RESULTS.MODULES.ASL = runtests('xASL_module_ASL_test');
% RESULTS.MODULES.POPULATION = runtests('xASL_module_Population_test');


%% RUN TESTS: SUBMODULES









end


