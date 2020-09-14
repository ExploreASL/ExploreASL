function [name,module,submodule,passed,tests] = ExploreASL_Unit_Test_Template
%ExploreASL_Unit_Test_Template Individual unit test template
%
% INPUT:        name   - Name of the module or submodule which is tested
%               module - Define if module or submodule is tested (true or false)
%
% OUTPUT:       test   - Structure with individual subtest results
%               passed - Result of all subtests combined (true or false)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  ...
%
% EXAMPLE:      [name,module,submodule,passed,tests] = ExploreASL_Unit_Test_Template;
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2020 ExploreASL

%% Initialize test structure (this does not have to be changed)
name = 'Template Test';
module = true;
submodule = false;

%% Test run 1
tests(1).testname = 'Check A';

% Start the test
testTime = tic;

% Run your test here
% ...

% Define one or multiple test conditions here
testCondition = true;

% Get test duration
tests(1).duration = toc(testTime);

% Evaluate your test
tests(1).passed = testCondition;




%% Test run 2
tests(2).testname = 'Check B';

% Start the test
testTime = tic;

% Run your test here
% ...

% Define one or multiple test conditions here
testCondition = true;

% Get test duration
tests(2).duration = toc(testTime);

% Evaluate your test
tests(2).passed = testCondition;





%% End of testing

% Check if an individual subtest failed
passed = true;
for it = 1:numel(tests)
    if ~tests(it).passed
        passed = false;
    end
end











