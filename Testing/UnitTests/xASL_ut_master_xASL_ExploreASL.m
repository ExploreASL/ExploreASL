function UnitTest = xASL_ut_master_xASL_ExploreASL(TestRepository)
%xASL_ut_master_xASL_ExploreASL Individual unit test for ExploreASL
%
% INPUT:        TestRepository - Path to test repository.
%
% OUTPUT:       UnitTest  - Test structure
%               name      - Name of tested module or submodule (char array)
%               unit      - Insert one of the following: 'Module', 'Submodule' or 'Function'
%               passed    - Result of all subtests combined (true or false)
%               test      - Structure with individual subtest results
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Should be run using xASL_ut_UnitTesting (!!!).
%
% EXAMPLE:      UnitTests(1) = xASL_ut_master_xASL_ExploreASL(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Initialize (from other directory)';

% Start the test
testTime = tic;

% Initialize test
currentDir = pwd;

% Make sure we're in the correct directory for testing
if ~strcmp(currentDir(end-9:end),'ExploreASL')
    error('Make sure to run the testing from xASL_test_UnitTesting and from the ExploreASL directory...');
end

% Initialize from correct directory
ExploreASL;

% Initialize again from other directory
cd('..');
x = ExploreASL;
cd(currentDir);

% Define one or multiple test conditions here
testCondition = true; % Fallback

% Check the basic fields first
if ~isstruct(x)
    testCondition = false;
end
if isfield(x, 'opts')
    if ~isfield(x.opts,'DatasetRoot'),     testCondition = false;      end
    if ~isfield(x.opts,'ImportModules'),   testCondition = false;      end
    if ~isfield(x.opts,'ProcessModules'),  testCondition = false;      end
    if ~isfield(x.opts,'bPause'),          testCondition = false;      end
    if ~isfield(x.opts,'iWorker'),         testCondition = false;      end
    if ~isfield(x.opts,'nWorkers'),        testCondition = false;      end
else
    testCondition = false;
end

% Now let's check the values
if isfield(x, 'opts')
    if isfield(x.opts,'DatasetRoot')
        if ~isempty(x.opts.DatasetRoot) || ~ischar(x.opts.DatasetRoot)
            testCondition = false;
        end
    end
    if isfield(x.opts,'ImportModules')
        if length(x.opts.ImportModules)<4 || sum(x.opts.ImportModules)>0 || ~isnumeric(x.opts.ImportModules)
            testCondition = false;
        end
    end
    if isfield(x.opts,'ProcessModules')
        if length(x.opts.ProcessModules)<3 || sum(x.opts.ProcessModules)>0 || ~isnumeric(x.opts.ProcessModules)
            testCondition = false;
        end
    end
    if isfield(x.opts,'bPause')
        if length(x.opts.bPause)>1 || sum(x.opts.bPause)>0 || ~isnumeric(x.opts.bPause)
            testCondition = false;
        end
    end
    if isfield(x.opts,'iWorker')
        if length(x.opts.iWorker)>1 || sum(x.opts.iWorker)>1 || ~isnumeric(x.opts.iWorker)
            testCondition = false;
        end
    end
    if isfield(x.opts,'nWorkers')
        if length(x.opts.nWorkers)>1 || sum(x.opts.nWorkers)>1 || ~isnumeric(x.opts.nWorkers)
            testCondition = false;
        end
    end
else
    testCondition = false;
end

% Clean-up
clearvars -except UnitTest TestRepository testCondition testTime

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


