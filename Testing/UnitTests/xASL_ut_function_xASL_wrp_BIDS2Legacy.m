function UnitTest = xASL_ut_function_xASL_wrp_BIDS2Legacy(TestRepository)
%xASL_ut_function_xASL_wrp_BIDS2Legacy Individual unit test for BIDS2Legacy
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_wrp_BIDS2Legacy(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright (c) 2015-2022 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Read in DRO test patient (default)';

% Start the test
testTime = tic;

% Define test patient paths
droTestPatientSource = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_3_0');
droTestPatient = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_3_0');
droSubject = 'sub-001'; % DRO subject
sessionNum = '_1';

% Copy test data to working directory
xASL_Copy(droTestPatientSource,droTestPatient);

% Initialize & run BIDS2Legacy
x = ExploreASL(droTestPatient, 0, 0, 0);

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~exist(fullfile(droTestPatient,'derivatives'),'dir')
    testCondition = false; % Test failed
end

% Check path
checkDirSubjSession = fullfile(droTestPatient,'derivatives','ExploreASL',[droSubject sessionNum]);

% Check ASL files
if ~exist(fullfile(checkDirSubjSession,'ASL_1','ASL4D.json'),'file') ...
    || ~exist(fullfile(checkDirSubjSession,'ASL_1','ASL4Dcontext_Source.tsv'),'file')
    testCondition = false; % Test failed
end
if ~xASL_exist(fullfile(checkDirSubjSession,'ASL_1','ASL4D.nii'),'file')
    testCondition = false; % Test failed
end

% Check M0 files
if ~exist(fullfile(checkDirSubjSession,'ASL_1','M0.json'),'file')
    testCondition = false; % Test failed
end
if ~xASL_exist(fullfile(checkDirSubjSession,'ASL_1','M0.nii')
    testCondition = false; % Test failed
end

% Check T1w files
if ~exist(fullfile(checkDirSubjSession,'T1.json'),'file')
    testCondition = false; % Test failed
end
if  ~xASL_exist(fullfile(checkDirSubjSession,'T1.nii')
    testCondition = false; % Test failed
end


% Remove the test files
xASL_delete(droTestPatient,true);

% Clean-Up
clearvars -except testTime TestRepository UnitTest testCondition

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
droTestPatientSource = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_3_0');
droTestPatient = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_3_0');

% Copy test data to working directory
xASL_Copy(droTestPatientSource,droTestPatient);

% The dataPar.json is not created for each subject anymore, which is why it
% was moved out of BIDS2Legacy. We need to create it manually here.
dataParJSON.x.dataset.subjectRegexp = '^sub-.*$';
dataParJSON.x.settings.Quality = 1;
xASL_adm_CreateDir(fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_3_0','derivatives','ExploreASL'));
fclose('all');
xASL_io_WriteJson(fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_3_0','derivatives','ExploreASL','dataPar.json'),dataParJSON);

% Initialize & run BIDS2Legacy
try
    x = ExploreASL(droTestPatient, 0, 0, 0);
catch
    x = false;
end

% Define one or multiple test conditions here
testCondition = true;

if ~exist(fullfile(droTestPatient,'derivatives'),'dir')
    testCondition = false; % Test failed
end
if ~isstruct(x)
    testCondition = false; % Test failed
end


% Remove the test files
xASL_delete(droTestPatient,true);

% Clean-Up
clearvars -except testTime TestRepository UnitTest testCondition

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);


end

