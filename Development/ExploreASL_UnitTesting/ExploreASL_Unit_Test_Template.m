function UnitTest = ExploreASL_Unit_Test_Template
%ExploreASL_Unit_Test_Template Individual unit test template
%
% INPUT:        n/a
%
% OUTPUT:       UnitTest  - Test structure
%               name      - Name of tested module or submodule (char array)
%               module    - True if module test
%               submodule - True if submodule test
%               passed    - Result of all subtests combined (true or false)
%               test      - Structure with individual subtest results
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  ...
%
% EXAMPLE:      [name,module,submodule,passed,tests] = ExploreASL_Unit_Test_Template;
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2020 ExploreASL

%% Initialize test structure (this does not have to be changed)
UnitTest.name = 'Template Test';
UnitTest.module = true;
UnitTest.submodule = false;

%% Test run 1
UnitTest.tests(1).testname = 'Check A';

% Start the test
testTime = tic;

% Run your test here
% ...

% Define one or multiple test conditions here
testCondition = true;

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;




%% Test run 2
UnitTest.tests(2).testname = 'Check B';

% Start the test
testTime = tic;

% Run your test here
% ...

% Define one or multiple test conditions here
testCondition = true;

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;





%% End of testing

% Check if an individual subtest failed
passed = true;
for it = 1:numel(UnitTest.tests)
    if ~UnitTest.tests(it).passed
        passed = false;
    end
end











