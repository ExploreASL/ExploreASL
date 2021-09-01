function UnitTest = xASL_ut_function_xASL_bids_BIDS2Legacy(TestRepository)
%xASL_ut_function_xASL_bids_BIDS2Legacy Individual unit test for BIDS2Legacy
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_bids_BIDS2Legacy(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Read in DRO test patient (default)';

% Start the test
testTime = tic;

% Define test patient paths
droTestPatientSource = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_2_0');
droTestPatient = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_2_0');
droSubject = 'sub-Sub1'; % DRO subject

% Copy test data to working directory
xASL_Copy(droTestPatientSource,fullfile(droTestPatient,'rawdata',droSubject),1);

% Prepare DRO
testVersion = '1.2.3';
xASL_bids_DRO2BIDS(droTestPatient,[],[],testVersion);

% Run BIDS2Legacy
xASL_bids_BIDS2Legacy(droTestPatient);

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~exist(fullfile(droTestPatient,'derivatives'),'dir')
    testCondition = false; % Test failed
end
if ~exist(fullfile(droTestPatient,'derivatives','ExploreASL','dataPar.json'),'file')
    testCondition = false; % Test failed
end

% Check ASL files
if ~exist(fullfile(droTestPatient,'derivatives','ExploreASL',droSubject,'ASL_1','ASL4D.json'),'file') ...
    || ~exist(fullfile(droTestPatient,'derivatives','ExploreASL',droSubject,'ASL_1','ASL4D_Source_aslcontext.tsv'),'file')
    testCondition = false; % Test failed
end
if ~exist(fullfile(droTestPatient,'derivatives','ExploreASL',droSubject,'ASL_1','ASL4D.nii'),'file') ...
    && ~exist(fullfile(droTestPatient,'derivatives','ExploreASL',droSubject,'ASL_1','ASL4D.nii.gz'),'file')
    testCondition = false; % Test failed
end

% Check M0 files
if ~exist(fullfile(droTestPatient,'derivatives','ExploreASL',droSubject,'ASL_1','M0.json'),'file')
    testCondition = false; % Test failed
end
if ~exist(fullfile(droTestPatient,'derivatives','ExploreASL',droSubject,'ASL_1','M0.nii'),'file') ...
    && ~exist(fullfile(droTestPatient,'derivatives','ExploreASL',droSubject,'ASL_1','M0.nii.gz'),'file')
    testCondition = false; % Test failed
end

% Check T1w files
if ~exist(fullfile(droTestPatient,'derivatives','ExploreASL',droSubject,'T1.json'),'file')
    testCondition = false; % Test failed
end
if  ~exist(fullfile(droTestPatient,'derivatives','ExploreASL',droSubject,'T1.nii'),'file') ...
    &&  ~exist(fullfile(droTestPatient,'derivatives','ExploreASL',droSubject,'T1.nii.gz'),'file')
    testCondition = false; % Test failed
end


% Remove the test files
xASL_delete(droTestPatient,true);

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% Test run 2

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Load DRO test patient in xASL (default)';

% Start the test
testTime = tic;

% Define test patient paths
droTestPatientSource = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_2_0');
droTestPatient = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_2_0');

% DRO subject
droSubject = 'sub-Sub1';

% Copy test data to working directory
xASL_Copy(droTestPatientSource,fullfile(droTestPatient,'rawdata',droSubject));

% Prepare DRO
testVersion = '1.2.3';
xASL_bids_DRO2BIDS(droTestPatient,[],[],testVersion);

% Run BIDS2Legacy
xASL_bids_BIDS2Legacy(droTestPatient);

% Initialize dataset
try
    [x] = ExploreASL_Master(droTestPatient, 0, 0);
catch
    x = false;
end

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~exist(fullfile(droTestPatient,'derivatives'),'dir')
    testCondition = false; % Test failed
end
if ~isstruct(x)
    testCondition = false; % Test failed
end


% Remove the test files
xASL_delete(droTestPatient,true);

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);


end

