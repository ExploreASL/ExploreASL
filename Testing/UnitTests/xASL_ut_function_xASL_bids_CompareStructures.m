function UnitTest = xASL_ut_function_xASL_bids_CompareStructures(TestRepository)
%xASL_ut_function_xASL_bids_CompareStructures Individual unit test for xASL_bids_CompareStructures
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_bids_CompareStructures(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Check for file and folder existence (identical)';

% Start the test
testTime = tic;

% Prepare the test
workingDir = fullfile(TestRepository,'UnitTesting','working_directory');
pathDatasetA = fullfile(workingDir,'datasetA');
pathDatasetB = fullfile(workingDir,'datasetB');
bPrintReport = [];
threshRmseNii = [];
detailedOutput = [];
printWarnings = [];
ignoreLogs = [];

% Create test directories
xASL_adm_CreateDir(pathDatasetA);
xASL_adm_CreateDir(pathDatasetB);

% Files & folders in A
xASL_adm_CreateDir(fullfile(pathDatasetA,'SubFolder1'));
xASL_adm_CreateDir(fullfile(pathDatasetA,'SubFolder2'));
fid1 = fopen(fullfile(workingDir,'datasetA','test_1.txt'), 'wt');
fid2 = fopen(fullfile(workingDir,'datasetA','test_2.txt'), 'wt');
fid3 = fopen(fullfile(workingDir,'datasetA','test_3.txt'), 'wt');
fclose('all');

% Files & folders in B
xASL_adm_CreateDir(fullfile(pathDatasetB,'SubFolder1'));
xASL_adm_CreateDir(fullfile(pathDatasetB,'SubFolder2'));
fid1 = fopen(fullfile(workingDir,'datasetB','test_1.txt'), 'wt');
fid2 = fopen(fullfile(workingDir,'datasetB','test_2.txt'), 'wt');
fid3 = fopen(fullfile(workingDir,'datasetB','test_3.txt'), 'wt');
fclose('all');

% Run your test here
[identical,results,reportTable] = xASL_bids_CompareStructures(...
    pathDatasetA,pathDatasetB,...
    bPrintReport,threshRmseNii,...
    detailedOutput,printWarnings,ignoreLogs);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~identical
    testCondition = false;
end
if ~isstruct(results)
    testCondition = false;
end
if ~istable(reportTable)
    testCondition = false;
end

% Clean-up after testing
xASL_delete(pathDatasetA,1);
xASL_delete(pathDatasetB,1);

% Clean-up
clearvars -except UnitTest TestRepository testCondition testTime

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;



%% Test run 2

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Check for file and folder existence (not identical)';

% Start the test
testTime = tic;

% Prepare the test
workingDir = fullfile(TestRepository,'UnitTesting','working_directory');
pathDatasetA = fullfile(workingDir,'datasetA');
pathDatasetB = fullfile(workingDir,'datasetB');
bPrintReport = [];
threshRmseNii = [];
detailedOutput = [];
printWarnings = [];
ignoreLogs = [];

% Create test directories
xASL_adm_CreateDir(pathDatasetA);
xASL_adm_CreateDir(pathDatasetB);

% Files & folders in A
xASL_adm_CreateDir(fullfile(pathDatasetA,'SubFolder1'));
xASL_adm_CreateDir(fullfile(pathDatasetA,'SubFolder2'));
xASL_adm_CreateDir(fullfile(pathDatasetA,'SubFolder1A'));
xASL_adm_CreateDir(fullfile(pathDatasetA,'SubFolder2A'));
fid1 = fopen(fullfile(workingDir,'datasetA','test_1.txt'), 'wt');
fid2 = fopen(fullfile(workingDir,'datasetA','test_2.txt'), 'wt');
fid3 = fopen(fullfile(workingDir,'datasetA','test_3.txt'), 'wt');
fid1 = fopen(fullfile(workingDir,'datasetA','test_1A.txt'), 'wt');
fid2 = fopen(fullfile(workingDir,'datasetA','test_2A.txt'), 'wt');
fid3 = fopen(fullfile(workingDir,'datasetA','test_3A.txt'), 'wt');
fclose('all');

% Files & folders in B
xASL_adm_CreateDir(fullfile(pathDatasetB,'SubFolder1'));
xASL_adm_CreateDir(fullfile(pathDatasetB,'SubFolder2'));
xASL_adm_CreateDir(fullfile(pathDatasetB,'SubFolder1B'));
xASL_adm_CreateDir(fullfile(pathDatasetB,'SubFolder2B'));
fid1 = fopen(fullfile(workingDir,'datasetB','test_1.txt'), 'wt');
fid2 = fopen(fullfile(workingDir,'datasetB','test_2.txt'), 'wt');
fid3 = fopen(fullfile(workingDir,'datasetB','test_3.txt'), 'wt');
fid1 = fopen(fullfile(workingDir,'datasetB','test_1B.txt'), 'wt');
fid2 = fopen(fullfile(workingDir,'datasetB','test_2B.txt'), 'wt');
fid3 = fopen(fullfile(workingDir,'datasetB','test_3B.txt'), 'wt');
fclose('all');

% Run your test here
[identical,results,reportTable] = xASL_bids_CompareStructures(...
    pathDatasetA,pathDatasetB,...
    bPrintReport,threshRmseNii,...
    detailedOutput,printWarnings,ignoreLogs);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if identical
    testCondition = false;
end
if ~isstruct(results)
    testCondition = false;
end
if ~istable(reportTable)
    testCondition = false;
end

% Clean-up after testing
xASL_delete(pathDatasetA,1);
xASL_delete(pathDatasetB,1);

% Clean-up
clearvars -except UnitTest TestRepository testCondition testTime

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;


%% Test run 3

% Give your individual subtest a name
UnitTest.tests(3).testname = 'Check JSON, TSV and general TEXT content (identical & not identical)';

% Start the test
testTime = tic;

% Prepare the test
workingDir = fullfile(TestRepository,'UnitTesting','working_directory');
pathDatasetA = fullfile(workingDir,'datasetA');
pathDatasetB = fullfile(workingDir,'datasetB');
bPrintReport = [];
threshRmseNii = [];
detailedOutput = 1;
printWarnings = 0;
ignoreLogs = 0;

% Test JSON A
testJSONa.fieldA = 'this is text field 1 in A';
testJSONa.fieldB = 'this is text field 2 in A';
testJSONa.fieldC = 'this is text field 3 in A';
testJSONa.fieldD = 123;
testJSONa.fieldE = [1, 2, 3];
testJSONa.fieldF = [1.23, 4.56, 7.89];

% Test JSON B
testJSONb.fieldA = 'this is text field 1 in B';
testJSONb.fieldB = 'this is text field 2 in B';
testJSONb.fieldC = 'this is text field 3 in B';
testJSONb.fieldD = 1234;
testJSONb.fieldE = [1, 2, 3, 4];
testJSONb.fieldF = [1.23, 4.56, 7.89, 10.11];

% Test TSV A
testTSVa = {'A', 'B', 'C'; 1, 2, 3; '[1, 2, 3]', '[4, 5, 6]', '[7, 8, 9]'};

% Test TSV B
testTSVb = {'D', 'E', 'F'; 1, 2, 3; '[1, 2, 3]', '[10, 11, 12]', '[7, 8, 9]'};

% Create test directories
xASL_adm_CreateDir(pathDatasetA);
xASL_adm_CreateDir(pathDatasetB);

% Files & folders in A
xASL_adm_CreateDir(fullfile(pathDatasetA,'SubFolder1'));
xASL_adm_CreateDir(fullfile(pathDatasetA,'SubFolder1A'));
fid1 = fopen(fullfile(workingDir,'datasetA','test_1.txt'), 'wt');
fprintf(fid1, 'Not really identical text in file 1');
fid2 = fopen(fullfile(workingDir,'datasetA','test_2.txt'), 'wt');
fprintf(fid2, 'Identical text in file 2');
fid1 = fopen(fullfile(workingDir,'datasetA','test_1A.txt'), 'wt');
fid2 = fopen(fullfile(workingDir,'datasetA','test_2A.txt'), 'wt');
fclose('all');

% Files & folders in B
xASL_adm_CreateDir(fullfile(pathDatasetB,'SubFolder1'));
xASL_adm_CreateDir(fullfile(pathDatasetB,'SubFolder1B'));
fid1 = fopen(fullfile(workingDir,'datasetB','test_1.txt'), 'wt');
fprintf(fid1, 'Not identical text in file 1');
fid2 = fopen(fullfile(workingDir,'datasetB','test_2.txt'), 'wt');
fprintf(fid2, 'Identical text in file 2');
fid1 = fopen(fullfile(workingDir,'datasetB','test_1B.txt'), 'wt');
fid2 = fopen(fullfile(workingDir,'datasetB','test_2B.txt'), 'wt');
fclose('all');

% Write identical JSON files
spm_jsonwrite(fullfile(workingDir,'datasetA','test_json_identical.txt'),testJSONa);
spm_jsonwrite(fullfile(workingDir,'datasetB','test_json_identical.txt'),testJSONa);

% Write not identical JSON files
spm_jsonwrite(fullfile(workingDir,'datasetA','test_json_not_identical.txt'),testJSONa);
spm_jsonwrite(fullfile(workingDir,'datasetB','test_json_not_identical.txt'),testJSONb);

% Write identical TSV files
xASL_tsvWrite(testTSVa,fullfile(workingDir,'datasetA','test_tsv_identical.tsv'));
xASL_tsvWrite(testTSVa,fullfile(workingDir,'datasetB','test_tsv_identical.tsv'));

% Write not identical TSV files
xASL_tsvWrite(testTSVa,fullfile(workingDir,'datasetA','test_tsv_not_identical.tsv'));
xASL_tsvWrite(testTSVb,fullfile(workingDir,'datasetB','test_tsv_not_identical.tsv'));

% Run your test here
[identical,results,reportTable] = xASL_bids_CompareStructures(...
    pathDatasetA,pathDatasetB,...
    bPrintReport,threshRmseNii,...
    detailedOutput,printWarnings,ignoreLogs);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if identical
    testCondition = false;
end
if ~isstruct(results)
    testCondition = false;
end
if ~istable(reportTable)
    testCondition = false;
end

% Check missing folders
if isfield(results,'datasetA') && isfield(results.datasetA,'missingFolders')
    if isempty(regexp(results.datasetA.missingFolders{1},'SubFolder1B', 'once'))
        testCondition = false;
    end
else
    testCondition = false;
end
if isfield(results,'datasetB') && isfield(results.datasetB,'missingFolders')
    if isempty(regexp(results.datasetB.missingFolders{1},'SubFolder1A', 'once'))
        testCondition = false;
    end
else
    testCondition = false;
end

% Check missing files
if isfield(results,'datasetA') && isfield(results.datasetA,'missingFiles')
    if isempty(regexp(results.datasetA.missingFiles{1},'test_1B.txt', 'once'))
        testCondition = false;
    end
    if isempty(regexp(results.datasetA.missingFiles{2},'test_2B.txt', 'once'))
        testCondition = false;
    end
else
    testCondition = false;
end
if isfield(results,'datasetB') && isfield(results.datasetB,'missingFiles')
    if isempty(regexp(results.datasetB.missingFiles{1},'test_1A.txt', 'once'))
        testCondition = false;
    end
    if isempty(regexp(results.datasetB.missingFiles{2},'test_2A.txt', 'once'))
        testCondition = false;
    end
else
    testCondition = false;
end

% Clean-up after testing
xASL_delete(pathDatasetA,1);
xASL_delete(pathDatasetB,1);

% Clean-up
clearvars -except UnitTest TestRepository testCondition testTime

% Get test duration
UnitTest.tests(3).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(3).passed = testCondition;


%% Test run 4

% Give your individual subtest a name
UnitTest.tests(4).testname = 'NIfTIs comparison (identical & not identical)';

% Start the test
testTime = tic;

% Prepare the test
workingDir = fullfile(TestRepository,'UnitTesting','working_directory');
pathDatasetA = fullfile(workingDir,'datasetA');
pathDatasetB = fullfile(workingDir,'datasetB');
bPrintReport = [];
threshRmseNii = [];
detailedOutput = 0;
printWarnings = 0;
ignoreLogs = 0;

% Create test directories
xASL_adm_CreateDir(pathDatasetA);
xASL_adm_CreateDir(pathDatasetB);

% Copy test NIfTIs
pathNIFTIperf = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_3_0_pre_release','sub-001','perf','sub-001_acq-001_asl.nii.gz');
pathNIFTIanat = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_3_0_pre_release','sub-001','anat','sub-001_acq-003_T1w.nii.gz');
xASL_Copy(pathNIFTIperf,fullfile(workingDir,'datasetA','testNifti.nii.gz'));
xASL_Copy(pathNIFTIperf,fullfile(workingDir,'datasetA','testNiftiIdentical.nii.gz'));
xASL_Copy(pathNIFTIanat,fullfile(workingDir,'datasetB','testNifti.nii.gz'));
xASL_Copy(pathNIFTIperf,fullfile(workingDir,'datasetB','testNiftiIdentical.nii.gz'));

% Run your test here
[identical,results,reportTable] = xASL_bids_CompareStructures(...
    pathDatasetA,pathDatasetB,...
    bPrintReport,threshRmseNii,...
    detailedOutput,printWarnings,ignoreLogs);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if identical
    testCondition = false;
end
if ~isstruct(results)
    testCondition = false;
end
if ~istable(reportTable)
    testCondition = false;
end
if size(results.differences,1)<1
    testCondition = false;
end

% Clean-up after testing
xASL_delete(pathDatasetA,1);
xASL_delete(pathDatasetB,1);

% Clean-up
clearvars -except UnitTest TestRepository testCondition testTime

% Get test duration
UnitTest.tests(4).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(4).passed = testCondition;



%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


