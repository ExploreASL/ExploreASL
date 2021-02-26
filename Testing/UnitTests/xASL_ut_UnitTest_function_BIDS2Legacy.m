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

% DRO subject
droSubject = 'sub-Sub1';

% Copy test data to working directory
xASL_Copy(droTestPatientSource,fullfile(droTestPatient,'rawdata',droSubject));

% Rename asl to perf
perfDirectory = fullfile(droTestPatient,'rawdata',droSubject,'perf');
xASL_Move(fullfile(droTestPatient,'rawdata',droSubject,'asl'),perfDirectory);
xASL_Move(fullfile(perfDirectory,'001_asl.json'),fullfile(perfDirectory,[droSubject,'_asl.json']));
xASL_Move(fullfile(perfDirectory,'001_asl.nii.gz'),fullfile(perfDirectory,[droSubject,'_asl.nii.gz']));
xASL_Move(fullfile(perfDirectory,'001_aslcontext.tsv'),fullfile(perfDirectory,[droSubject,'_aslcontext.tsv']));

% Run BIDS2Legacy

% Create dummy dataset_description.json
json = struct;
json.Name = "DRO_Digital_Reference_Object";
json.BIDSVersion = "1.5.0";
json.DatasetType = "raw";
json.License = "RandomText";
json.Authors = "RandomText";
json.Acknowledgements = "RandomText";
json.HowToAcknowledge = "Please cite this paper: https://www.ncbi.nlm.nih.gov/pubmed/001012092119281";
json.Funding = "RandomText";
json.EthicsApprovals = "RandomText";
json.ReferencesAndLinks = "RandomText";
json.DatasetDOI = "RandomText";

% Write file
spm_jsonwrite(fullfile(droTestPatient,'rawdata','dataset_description.json'),json);

% Run BIDS2Legacy
xASL_bids_BIDS2Legacy(droTestPatient);

% Define one or multiple test conditions here
testCondition = true; % Fallback
if ~exist(fullfile(droTestPatient,'derivatives'),'dir')
    testCondition = false; % Test failed
end
if ~exist(fullfile(droTestPatient,'derivatives','ExploreASL','DataPar.json'),'file')
    testCondition = false; % Test failed
end
if ~exist(fullfile(droTestPatient,'derivatives','ExploreASL',droSubject,'ASL_1','ASL4D.json'),'file')
    testCondition = false; % Test failed
end
if ~exist(fullfile(droTestPatient,'derivatives','ExploreASL',droSubject,'ASL_1','ASL4D.nii.gz'),'file')
    testCondition = false; % Test failed
end
if ~exist(fullfile(droTestPatient,'derivatives','ExploreASL',droSubject,'ASL_1','ASL4D_aslcontext.tsv'),'file')
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

