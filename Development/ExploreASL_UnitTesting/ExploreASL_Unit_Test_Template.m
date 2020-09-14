function [test,passed] = ExploreASL_Unit_Test_Template(name,module)
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
% EXAMPLE:      [test,passed] = ExploreASL_Unit_Test_Template('Template Test',true);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2020 ExploreASL

%% Initialize test structure (this does not have to be changed)
test.name = [];
test.module = [];
test.submodule = [];
test.testname = [];

%% Test run 1
test(1).testname = 'Check A';

% Start the test
testTime = tic;

% Run your test here
% ...

% Define one or multiple test conditions here
testCondition = true;

% Get test duration
test(1).duration = toc(testTime);

% Evaluate your test
test(1).passed = testCondition;




%% Test run 2
test(2).testname = 'Check B';

% Start the test
testTime = tic;

% Run your test here
% ...

% Define one or multiple test conditions here
testCondition = true;

% Get test duration
test(2).duration = toc(testTime);

% Evaluate your test
test(2).passed = testCondition;





%% End of testing

% Assignment of previously defined fields
passed = true;
for it = 1:numel(test)
    test(it).name = name;
    test(it).module = module;
    test(it).submodule = ~module;
    % Check if an individual subtest failed
    if ~test(it).passed
        passed = false;
    end
end











