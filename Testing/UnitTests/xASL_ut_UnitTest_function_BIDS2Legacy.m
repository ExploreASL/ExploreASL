function UnitTest = xASL_ut_UnitTest_function_BIDS2Legacy(TestRepository)
%xASL_ut_UnitTest_function_BIDS2Legacy Individual unit test for BIDS2Legacy
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
% EXAMPLE:      UnitTests(1) = xASL_ut_UnitTest_function_BIDS2Legacy(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL

%% Initialize test structure

% Insert test name here
UnitTest.name = 'xASL_bids_BIDS2Legacy';

% Define whether you are testing a module, submodule or function
UnitTest.unit = 'Function';

%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Read in DRO test patient (default)';

% Start the test
testTime = tic;

% Define test patient paths
droTestPatientSource = fullfile(TestRepository,'UnitTesting','dro_files','test_patient');
droTestPatient = fullfile(TestRepository,'UnitTesting','working_directory','test_patient');

% Copy test data to working directory
copyfile(droTestPatientSource,fullfile(droTestPatient,'rawdata'));
    
% Run BIDS2Legacy
warning('We need to add the dataset_description.json to the rawdata folder first...');
% xASL_bids_BIDS2Legacy(droTestPatient);

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~exist(fullfile(droTestPatient,'sourcedata'),'dir')
    testCondition = false; % Test failed
end

% Remove the test files
xASL_delete(droTestPatient,true);

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;

%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);


end

