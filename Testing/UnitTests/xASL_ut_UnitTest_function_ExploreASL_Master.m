function UnitTest = xASL_ut_UnitTest_function_ExploreASL_Master(TestRepository)
%xASL_ut_UnitTest_function_ExploreASL_Master Individual unit test for
%ExploreASL_Master
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
% DESCRIPTION:  Should be run using xASL_ut_UnitTesting.
%
% EXAMPLE:      UnitTests(1) = xASL_ut_UnitTest_function_ExploreASL_Master(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL

%% Initialize test structure

% Insert test name here
UnitTest.name = 'ExploreASL_Master';

% Define whether you are testing a module, submodule or function
UnitTest.unit = 'Function';

%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Initialize (without arguments)';

% Start the test
testTime = tic;

% Read test files
[x] = ExploreASL_Master();

% Define one or multiple test conditions here
testCondition = true; % Fallback

% Check the basic fields first
if ~isstruct(x)
    testCondition = false;
end
if ~isfield(x,'DataParPath') 
    testCondition = false;
end
if ~isfield(x,'ImportArray')
    testCondition = false;
end
if ~isfield(x,'ProcessArray')
    testCondition = false;
end
if ~isfield(x,'SkipPause')
    testCondition = false;
end
if ~isfield(x,'iWorker')
    testCondition = false;
end
if ~isfield(x,'nWorkers')
    testCondition = false;
end

% Now let's check the values
if isfield(x,'DataParPath')
    if ~isempty(x.DataParPath) || ~ischar(x.DataParPath)
        testCondition = false;
    end
end
if isfield(x,'ImportArray')
    if length(x.ImportArray)<4 || sum(x.ImportArray)>0 || ~isnumeric(x.ImportArray)
        testCondition = false;
    end
end
if isfield(x,'ProcessArray')
    if length(x.ProcessArray)<3 || sum(x.ProcessArray)>0 || ~isnumeric(x.ProcessArray)
        testCondition = false;
    end
end
if isfield(x,'SkipPause')
    if length(x.SkipPause)>1 || sum(x.SkipPause)>0 || ~isnumeric(x.SkipPause)
        testCondition = false;
    end
end
if isfield(x,'iWorker')
    if length(x.iWorker)>1 || sum(x.iWorker)>1 || ~isnumeric(x.iWorker)
        testCondition = false;
    end
end
if isfield(x,'nWorkers')
    if length(x.nWorkers)>1 || sum(x.nWorkers)>1 || ~isnumeric(x.nWorkers)
        testCondition = false;
    end
end

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


