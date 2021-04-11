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
if ~isfield(x,'DataParPath'),     testCondition = false;      end
if ~isfield(x,'ImportModules'),   testCondition = false;      end
if ~isfield(x,'ProcessModules'),  testCondition = false;      end
if ~isfield(x,'bPause'),          testCondition = false;      end
if ~isfield(x,'iWorker'),         testCondition = false;      end
if ~isfield(x,'nWorkers'),        testCondition = false;      end

% Now let's check the values
if isfield(x,'DataParPath')
    if ~isempty(x.DataParPath) || ~ischar(x.DataParPath)
        testCondition = false;
    end
end
if isfield(x,'ImportModules')
    if length(x.ImportModules)<4 || sum(x.ImportModules)>0 || ~isnumeric(x.ImportModules)
        testCondition = false;
    end
end
if isfield(x,'ProcessModules')
    if length(x.ProcessModules)<3 || sum(x.ProcessModules)>0 || ~isnumeric(x.ProcessModules)
        testCondition = false;
    end
end
if isfield(x,'bPause')
    if length(x.bPause)>1 || sum(x.bPause)>0 || ~isnumeric(x.bPause)
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


%% Test run 2

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Initialize (with arguments)';

% Start the test
testTime = tic;

% Read test files
[x] = ExploreASL_Master('',0,0,0,1,1);

% Define one or multiple test conditions here
testCondition = true; % Fallback

% Check the basic fields first
if ~isstruct(x)
    testCondition = false;
end
if ~isfield(x,'DataParPath'),     testCondition = false;      end
if ~isfield(x,'ImportModules'),   testCondition = false;      end
if ~isfield(x,'ProcessModules'),  testCondition = false;      end
if ~isfield(x,'bPause'),          testCondition = false;      end
if ~isfield(x,'iWorker'),         testCondition = false;      end
if ~isfield(x,'nWorkers'),        testCondition = false;      end

% Now let's check the values
if isfield(x,'DataParPath')
    if ~isempty(x.DataParPath) || ~ischar(x.DataParPath)
        testCondition = false;
    end
end
if isfield(x,'ImportModules')
    if length(x.ImportModules)<4 || sum(x.ImportModules)>0 || ~isnumeric(x.ImportModules)
        testCondition = false;
    end
end
if isfield(x,'ProcessModules')
    if length(x.ProcessModules)<3 || sum(x.ProcessModules)>0 || ~isnumeric(x.ProcessModules)
        testCondition = false;
    end
end
if isfield(x,'bPause')
    if length(x.bPause)>1 || x.bPause~=0 || ~isnumeric(x.bPause)
        testCondition = false;
    end
end
if isfield(x,'iWorker')
    if length(x.iWorker)>1 || x.iWorker~=1 || ~isnumeric(x.iWorker)
        testCondition = false;
    end
end
if isfield(x,'nWorkers')
    if length(x.nWorkers)>1 || x.nWorkers~=1 || ~isnumeric(x.nWorkers)
        testCondition = false;
    end
end

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;


%% Test run 3

% Give your individual subtest a name
UnitTest.tests(3).testname = 'Initialize (with arrays)';

% Start the test
testTime = tic;

% Read test files
[x] = ExploreASL_Master('',[0 0 0 0],[0 0 0 0],0,1,1);

% Define one or multiple test conditions here
testCondition = true; % Fallback

% Check the basic fields first
if ~isstruct(x)
    testCondition = false;
end
if ~isfield(x,'DataParPath'),     testCondition = false;      end
if ~isfield(x,'ImportModules'),   testCondition = false;      end
if ~isfield(x,'ProcessModules'),  testCondition = false;      end
if ~isfield(x,'bPause'),          testCondition = false;      end
if ~isfield(x,'iWorker'),         testCondition = false;      end
if ~isfield(x,'nWorkers'),        testCondition = false;      end

% Now let's check the values
if isfield(x,'DataParPath')
    if ~isempty(x.DataParPath) || ~ischar(x.DataParPath)
        testCondition = false;
    end
end
if isfield(x,'ImportModules')
    if length(x.ImportModules)<4 || sum(x.ImportModules)>0 || ~isnumeric(x.ImportModules)
        testCondition = false;
    end
end
if isfield(x,'ProcessModules')
    if length(x.ProcessModules)<3 || sum(x.ProcessModules)>0 || ~isnumeric(x.ProcessModules)
        testCondition = false;
    end
end
if isfield(x,'bPause')
    if length(x.bPause)>1 || x.bPause~=0 || ~isnumeric(x.bPause)
        testCondition = false;
    end
end
if isfield(x,'iWorker')
    if length(x.iWorker)>1 || x.iWorker~=1 || ~isnumeric(x.iWorker)
        testCondition = false;
    end
end
if isfield(x,'nWorkers')
    if length(x.nWorkers)>1 || x.nWorkers~=1 || ~isnumeric(x.nWorkers)
        testCondition = false;
    end
end

% Get test duration
UnitTest.tests(3).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(3).passed = testCondition;


%% Test run 4

% Give your individual subtest a name
UnitTest.tests(4).testname = 'Initialize (DCM2NIFTI)';

% Start the test
testTime = tic;

% Copy test patient
testPatientSource = fullfile(TestRepository,'UnitTesting','synthetic_dcm','test_patient'); 
testPatientDestination = fullfile(TestRepository,'UnitTesting','working_directory','test_patient');
xASL_Copy(testPatientSource, testPatientDestination);

% Read test files
[x] = ExploreASL_Master(fullfile(testPatientDestination,'sourceStructure.json'),[1 0 0 0],0,0,1,1);

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~exist(fullfile(testPatientDestination,'analysis'),'dir'),                   testCondition = false;          end
if ~exist(fullfile(testPatientDestination,'analysis','Sub1'),'dir'),            testCondition = false;          end
if ~exist(fullfile(testPatientDestination,'analysis','Sub1','ASL_1'),'dir'),    testCondition = false;          end
if ~exist(fullfile(testPatientDestination,'analysis','Sub1','T1w_1'),'dir'),    testCondition = false;          end
if ~exist(fullfile(testPatientDestination,'analysis','Sub1','ASL_1','ASL4D.json'),'file'),      testCondition = false;          end
if ~exist(fullfile(testPatientDestination,'analysis','Sub1','ASL_1','ASL4D.nii'),'file'),       testCondition = false;          end
if ~exist(fullfile(testPatientDestination,'analysis','Sub1','ASL_1','M0.json'),'file'),         testCondition = false;          end
if ~exist(fullfile(testPatientDestination,'analysis','Sub1','ASL_1','M0.nii'),'file'),          testCondition = false;          end
if ~exist(fullfile(testPatientDestination,'analysis','Sub1','T1w_1','T1w.json'),'file'),        testCondition = false;          end
if ~exist(fullfile(testPatientDestination,'analysis','Sub1','T1w_1','T1w.nii'),'file'),         testCondition = false;          end

% Delete test data
xASL_delete(testPatientDestination,true)

% Get test duration
UnitTest.tests(4).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(4).passed = testCondition;


%% Test run 5

% ...


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


