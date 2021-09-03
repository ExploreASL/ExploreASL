function UnitTest = xASL_ut_master_xASL_ExploreASL_Master(TestRepository)
%xASL_ut_master_xASL_ExploreASL_Master Individual unit test for ExploreASL_Master
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
% 1.  Initialize  (Without arguments)
% 2.  Initialize  (With arguments)
% 3.  Initialize  (With arrays)
% 4.  DRO 2.2.0   (DCM2NIFTI)
% 5.  DRO 2.2.0   (NII2BIDS)
% 6.  DRO 2.2.0   (Deface, BIDS2Legacy)
% 7.  DRO 2.2.0   (Deface, BIDS2Legacy with dataPar.json)
% 8.  DRO 2.2.0   (Run processing starting from derivatives with directory input)
% 9.  DRO 2.2.0   (Run processing starting from derivatives with dataPar.json input (outdated))
% 10. DRO 2.2.0   (Full pipeline, rawdata->results)
% 11. DRO 2.3.0   (Pre-release version, multi-session BIDS to legacy)
%
% EXAMPLE:      UnitTests(1) = xASL_ut_master_xASL_ExploreASL_Master(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Initialize (Without arguments)';

% Start the test
testTime = tic;

% Read test files
[x] = ExploreASL_Master();
testVersion = x.Version;

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
clearvars -except UnitTest TestRepository testCondition testTime testVersion

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% Test run 2

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Initialize (With arguments)';

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
clearvars -except UnitTest TestRepository testCondition testTime testVersion

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;


%% Test run 3

% Give your individual subtest a name
UnitTest.tests(3).testname = 'Initialize (With arrays)';

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
clearvars -except UnitTest TestRepository testCondition testTime testVersion

% Get test duration
UnitTest.tests(3).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(3).passed = testCondition;


%% Test run 4

% Give your individual subtest a name
UnitTest.tests(4).testname = 'DRO 2.2.0 (DCM2NIFTI)';

% Start the test
testTime = tic;

% Copy test patient
testPatientSource = fullfile(TestRepository,'UnitTesting','synthetic_dcm','test_patient_2_2_0'); 
testPatientDestination = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_2_0');
xASL_Copy(testPatientSource, testPatientDestination, 1);

% Fallback
testCondition = true;

% Read test files
try
    [x] = ExploreASL_Master(testPatientDestination,[1 0 0 0],0,0,1,1);
catch ME
    warning(ME.identifier, '%s', ME.message);
    testCondition = false;
    diary off;
end

% Define one or multiple test conditions here
if ~exist('x','var')
    testCondition = false;
else
    if ~isstruct(x)
        testCondition = false;
    end
end
if ~exist(fullfile(testPatientDestination,'temp'),'dir'),                   testCondition = false;          end
if ~exist(fullfile(testPatientDestination,'temp','Sub1'),'dir'),            testCondition = false;          end
if ~exist(fullfile(testPatientDestination,'temp','Sub1','ASL_1'),'dir'),    testCondition = false;          end
if ~exist(fullfile(testPatientDestination,'temp','Sub1','T1w_1'),'dir'),    testCondition = false;          end
if ~exist(fullfile(testPatientDestination,'temp','Sub1','ASL_1','ASL4D.json'),'file'),      testCondition = false;          end
if ~exist(fullfile(testPatientDestination,'temp','Sub1','ASL_1','ASL4D.nii'),'file'),       testCondition = false;          end
if ~exist(fullfile(testPatientDestination,'temp','Sub1','ASL_1','M0.json'),'file'),         testCondition = false;          end
if ~exist(fullfile(testPatientDestination,'temp','Sub1','ASL_1','M0.nii'),'file'),          testCondition = false;          end
if ~exist(fullfile(testPatientDestination,'temp','Sub1','T1w_1','T1w.json'),'file'),        testCondition = false;          end
if ~exist(fullfile(testPatientDestination,'temp','Sub1','T1w_1','T1w.nii'),'file'),         testCondition = false;          end

% Delete test data
xASL_delete(testPatientDestination,true)

% Clean-up
clearvars -except UnitTest TestRepository testCondition testTime testVersion

% Get test duration
UnitTest.tests(4).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(4).passed = testCondition;


%% Test run 5

% Give your individual subtest a name
UnitTest.tests(5).testname = 'DRO 2.2.0 (NII2BIDS)';

% Start the test
testTime = tic;

% Set-up DRO
droTestPatientSource = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_2_0');
droTestPatient = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_2_0');
testPatientDestination = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_2_0');
droSubject = 'sub-Sub1'; % DRO subject
xASL_Copy(droTestPatientSource,fullfile(droTestPatient,'rawdata',droSubject),1);
xASL_bids_DRO2BIDS(droTestPatient,[],[],testVersion); % Prepare DRO

% Convert BIDS back to temp for the testing
xASL_Move(fullfile(droTestPatient,'rawdata'),fullfile(droTestPatient,'temp'))
xASL_Move(fullfile(droTestPatient,'temp',droSubject),fullfile(droTestPatient,'temp','Sub1'))
xASL_Move(fullfile(droTestPatient,'temp','Sub1','anat'),fullfile(droTestPatient,'temp','Sub1','T1w_1'))
xASL_Move(fullfile(droTestPatient,'temp','Sub1','perf'),fullfile(droTestPatient,'temp','Sub1','ASL_1'))
xASL_delete(fullfile(droTestPatient,'temp','dataset_description.json'))
xASL_Move(fullfile(droTestPatient,'temp','Sub1','T1w_1','sub-Sub1_T1w.json'),fullfile(droTestPatient,'temp','Sub1','T1w_1','T1w.json'))
xASL_Move(fullfile(droTestPatient,'temp','Sub1','T1w_1','sub-Sub1_T1w.nii.gz'),fullfile(droTestPatient,'temp','Sub1','T1w_1','T1w.nii.gz'))
xASL_delete(fullfile(droTestPatient,'temp','Sub1','ASL_1','sub-Sub1_aslcontext.tsv'))
xASL_Move(fullfile(droTestPatient,'temp','Sub1','ASL_1','sub-Sub1_asl.json'),fullfile(droTestPatient,'temp','Sub1','ASL_1','ASL4D.json'))
xASL_Move(fullfile(droTestPatient,'temp','Sub1','ASL_1','sub-Sub1_asl.nii.gz'),fullfile(droTestPatient,'temp','Sub1','ASL_1','ASL4D.nii.gz'))

% Fallback
testCondition = true;

% Read test files
try
    [x] = ExploreASL_Master(testPatientDestination,[0 1 0 0],0,0,1,1);
catch ME
    warning(ME.identifier, '%s', ME.message);
    testCondition = false;
    diary off;
end

% Define one or multiple test conditions here
if ~exist('x','var')
    testCondition = false;
else
    if ~isstruct(x)
        testCondition = false;
    end
end

% Check ASL files
if ~exist(fullfile(droTestPatient,'rawdata',droSubject,'perf','sub-Sub1_asl.json'),'file') ...
    || ~exist(fullfile(droTestPatient,'rawdata',droSubject,'perf','sub-Sub1_aslcontext.tsv'),'file')
    testCondition = false; % Test failed
end
if ~exist(fullfile(droTestPatient,'rawdata',droSubject,'perf','sub-Sub1_asl.nii'),'file') ...
    && ~exist(fullfile(droTestPatient,'rawdata',droSubject,'perf','sub-Sub1_asl.nii.gz'),'file')
    testCondition = false; % Test failed
end

% Check T1w files
if ~exist(fullfile(droTestPatient,'rawdata',droSubject,'anat','sub-Sub1_T1w.json'),'file')
    testCondition = false; % Test failed
end
if ~exist(fullfile(droTestPatient,'rawdata',droSubject,'anat','sub-Sub1_T1w.nii'),'file') ...
    && ~exist(fullfile(droTestPatient,'rawdata',droSubject,'anat','sub-Sub1_T1w.nii.gz'),'file')
    testCondition = false; % Test failed
end

% Delete test data
xASL_delete(testPatientDestination,true)

% Clean-up
clearvars -except UnitTest TestRepository testCondition testTime testVersion

% Get test duration
UnitTest.tests(5).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(5).passed = testCondition;


%% Test run 6

% Give your individual subtest a name
UnitTest.tests(6).testname = 'DRO 2.2.0 (Deface, BIDS2Legacy)';

% Start the test
testTime = tic;

% Set-up DRO
droTestPatientSource = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_2_0');
droTestPatient = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_2_0');
testPatientDestination = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_2_0');
droSubject = 'sub-Sub1'; % DRO subject
xASL_Copy(droTestPatientSource,fullfile(droTestPatient,'rawdata',droSubject),1);
xASL_bids_DRO2BIDS(droTestPatient,[],[],testVersion); % Prepare DRO

% Fallback
testCondition = true;

% Read test files
try
    [x] = ExploreASL_Master(testPatientDestination,[0 0 1 1],0,0,1,1);
catch ME
    warning(ME.identifier, '%s', ME.message);
    testCondition = false;
    diary off;
end

% Define one or multiple test conditions here
if ~exist('x','var')
    testCondition = false;
else
    if ~isstruct(x)
        testCondition = false;
    end
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

% Check T1w files
if ~exist(fullfile(droTestPatient,'derivatives','ExploreASL',droSubject,'T1.json'),'file')
    testCondition = false; % Test failed
end
if  ~exist(fullfile(droTestPatient,'derivatives','ExploreASL',droSubject,'T1.nii'),'file') ...
    &&  ~exist(fullfile(droTestPatient,'derivatives','ExploreASL',droSubject,'T1.nii.gz'),'file')
    testCondition = false; % Test failed
end

% Delete test data
xASL_delete(testPatientDestination,true)

% Clean-up
clearvars -except UnitTest TestRepository testCondition testTime testVersion

% Get test duration
UnitTest.tests(6).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(6).passed = testCondition;


%% Test run 7

% Give your individual subtest a name
UnitTest.tests(7).testname = 'DRO 2.2.0 (Deface, BIDS2Legacy with dataPar.json)';

% Start the test
testTime = tic;

% Set-up DRO
droTestPatientSource = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_2_0');
droTestPatient = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_2_0');
testPatientDestination = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_2_0');
droSubject = 'sub-Sub1'; % DRO subject
xASL_Copy(droTestPatientSource,fullfile(droTestPatient,'rawdata',droSubject),1);
xASL_bids_DRO2BIDS(droTestPatient,[],[],testVersion); % Prepare DRO
% Create dataPar.json
dataParStruct.x.settings.Quality = 0;
dataParStruct.x.S.Atlases = {'TotalGM','DeepWM','Hammers','HOcort_CONN','HOsub_CONN','Mindboggle_OASIS_DKT31_CMA'};
spm_jsonwrite(fullfile(droTestPatient,'dataPar.json'),dataParStruct);

% Fallback
testCondition = true;

% Read test files
try
    [x] = ExploreASL_Master(testPatientDestination,[0 0 1 1],0,0,1,1);
catch ME
    warning(ME.identifier, '%s', ME.message);
    testCondition = false;
    diary off;
end

% Define one or multiple test conditions here
if ~exist('x','var')
    testCondition = false;
else
    if ~isstruct(x)
        testCondition = false;
    end
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

% Check T1w files
if ~exist(fullfile(droTestPatient,'derivatives','ExploreASL',droSubject,'T1.json'),'file')
    testCondition = false; % Test failed
end
if  ~exist(fullfile(droTestPatient,'derivatives','ExploreASL',droSubject,'T1.nii'),'file') ...
    &&  ~exist(fullfile(droTestPatient,'derivatives','ExploreASL',droSubject,'T1.nii.gz'),'file')
    testCondition = false; % Test failed
end

% Check dataPar.json
if exist(fullfile(droTestPatient,'derivatives','ExploreASL','dataPar.json'),'file')
    testContent = spm_jsonread(fullfile(droTestPatient,'derivatives','ExploreASL','dataPar.json'));
    if isfield(testContent,'x')
        if isfield(testContent.x,'S')
            if isfield(testContent.x.S,'Atlases')
                if ~(numel(testContent.x.S.Atlases)==6)
                    testCondition = false; % Test failed
                end
            else
                testCondition = false; % Test failed
            end
        else
            testCondition = false; % Test failed
        end
    else
        testCondition = false; % Test failed
    end
else
    testCondition = false; % Test failed
end

% Delete test data
xASL_delete(testPatientDestination,true)

% Clean-up
clearvars -except UnitTest TestRepository testCondition testTime testVersion

% Get test duration
UnitTest.tests(7).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(7).passed = testCondition;


%% Test run 8

% Give your individual subtest a name
UnitTest.tests(8).testname = 'DRO 2.2.0 (Run processing starting from derivatives with directory input)';

% Start the test
testTime = tic;

% Set-up DRO
droTestPatientSource = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_2_0');
droTestPatient = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_2_0');
testPatientDestination = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_2_0');
droSubject = 'sub-Sub1'; % DRO subject
xASL_Copy(droTestPatientSource,fullfile(droTestPatient,'rawdata',droSubject),1);
xASL_bids_DRO2BIDS(droTestPatient,[],[],testVersion); % Prepare DRO
% Create dataPar.json
dataParStruct.x.settings.Quality = 0;
dataParStruct.x.S.Atlases = {'TotalGM','DeepWM','Hammers','HOcort_CONN','HOsub_CONN','Mindboggle_OASIS_DKT31_CMA'};
spm_jsonwrite(fullfile(droTestPatient,'dataPar.json'),dataParStruct);

% Fallback
testCondition = true;

% Prepare the derivatives data
try
    [x] = ExploreASL_Master(testPatientDestination,[0 0 1 1],0,0,1,1);
catch ME
    warning(ME.identifier, '%s', ME.message);
    testCondition = false;
    diary off;
end

% Actual test: run processing starting from derivatives with directory input
try
    [x] = ExploreASL_Master(testPatientDestination,0,1,0,1,1);
catch ME
    warning(ME.identifier, '%s', ME.message);
    testCondition = false;
    diary off;
end

% Test conditions
if ~exist('x','var')
    testCondition = false;
else
    if ~isstruct(x)
        testCondition = false;
    end
end

% Check directories
testDirsAndFiles.derivativesDir = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_2_0','derivatives');
testDirsAndFiles.exploreASLDir = fullfile(testDirsAndFiles.derivativesDir,'ExploreASL');
testDirsAndFiles.populationDir = fullfile(testDirsAndFiles.exploreASLDir,'Population');
testDirsAndFiles.subDir = fullfile(testDirsAndFiles.exploreASLDir,'sub-Sub1');
testDirsAndFiles.aslDir = fullfile(testDirsAndFiles.subDir,'ASL_1');

% Check files
testDirsAndFiles.catReport = fullfile(testDirsAndFiles.subDir,'catreport_T1.pdf');
testDirsAndFiles.aslReport = fullfile(testDirsAndFiles.subDir,'xASL_Report_sub-Sub1.pdf');

% Iterate over test directories and files
fieldsTestDirsAndFiles = fieldnames(testDirsAndFiles);
for iField = 1:numel(fieldsTestDirsAndFiles)
    if ~xASL_exist(testDirsAndFiles.(fieldsTestDirsAndFiles{iField}),'file') && ~xASL_exist(testDirsAndFiles.(fieldsTestDirsAndFiles{iField}),'dir')
        testCondition = false;
    end
end

% Delete test data
xASL_delete(testPatientDestination,true)

% Clean-up
clearvars -except UnitTest TestRepository testCondition testTime testVersion

% Get test duration
UnitTest.tests(8).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(8).passed = testCondition;


%% Test run 9

% Give your individual subtest a name
UnitTest.tests(9).testname = 'DRO 2.2.0 (Run processing starting from derivatives with dataPar.json input (outdated))';

% Start the test
testTime = tic;

% Set-up DRO
droTestPatientSource = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_2_0');
droTestPatient = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_2_0');
testPatientDestination = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_2_0');
droSubject = 'sub-Sub1'; % DRO subject
xASL_Copy(droTestPatientSource,fullfile(droTestPatient,'rawdata',droSubject),1);
xASL_bids_DRO2BIDS(droTestPatient,[],[],testVersion); % Prepare DRO
% Create dataPar.json
dataParStruct.x.settings.Quality = 0;
dataParStruct.x.S.Atlases = {'TotalGM','DeepWM','Mindboggle_OASIS_DKT31_CMA'};
spm_jsonwrite(fullfile(droTestPatient,'dataPar.json'),dataParStruct);

% Fallback
testCondition = true;

% Prepare the derivatives data
try
    [x] = ExploreASL_Master(testPatientDestination,[0 0 1 1],0,0,1,1);
catch ME
    warning(ME.identifier, '%s', ME.message);
    testCondition = false;
    diary off;
end

% Actual test: run processing starting from derivatives with dataPar.json input (outdated)
try
    [x] = ExploreASL_Master(fullfile(testPatientDestination,'derivatives','ExploreASL','dataPar.json'),0,1,0,1,1);
catch ME
    warning(ME.identifier, '%s', ME.message);
    testCondition = false;
    diary off;
end

% Add test conditions here ...
if ~exist('x','var')
    testCondition = false;
else
    if ~isstruct(x)
        testCondition = false;
    end
end

% Delete test data
xASL_delete(testPatientDestination,true)

% Clean-up
clearvars -except UnitTest TestRepository testCondition testTime testVersion

% Get test duration
UnitTest.tests(9).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(9).passed = testCondition;


%% Test run 10

% Give your individual subtest a name
UnitTest.tests(10).testname = 'DRO 2.2.0 (Full pipeline, rawdata->results)';

% Start the test
testTime = tic;

% Set-up DRO
droTestPatientSource = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_2_0');
droTestPatient = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_2_0');
testPatientDestination = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_2_0');
droSubject = 'sub-Sub1'; % DRO subject
xASL_Copy(droTestPatientSource,fullfile(droTestPatient,'rawdata',droSubject),1);
xASL_bids_DRO2BIDS(droTestPatient,droSubject,false,testVersion); % Prepare DRO, keep ground truth data for comparison

% Fallback
testCondition = true;

% Read test files
try
    [x] = ExploreASL_Master(testPatientDestination,[0 0 1 1],1,0,1,1);
catch ME
    warning(ME.identifier, '%s', ME.message);
    testCondition = false;
    diary off;
end

% Define one or multiple test conditions here
if ~exist('x','var')
    testCondition = false;
else
    if ~isstruct(x)
        testCondition = false;
    end
end

% Check files and folders
if ~exist(fullfile(droTestPatient,'derivatives','ExploreASL',droSubject),'dir')
    testCondition = false; % Test failed
end
if ~exist(fullfile(droTestPatient,'derivatives','ExploreASL','dataPar.json'),'file')
    testCondition = false; % Test failed
end

% ...

% Compare image data
groundTruthM0File = fullfile(droTestPatient,'rawdata',droSubject,'ground_truth','003_ground_truth_m0.nii');
groundTruthPerfusionFile = fullfile(droTestPatient,'rawdata',droSubject,'ground_truth','003_ground_truth_perfusion_rate.nii');
derivedM0File = fullfile(droTestPatient,'derivatives','ExploreASL',droSubject,'ASL_1','M0.nii');
derivedCBFFile = fullfile(droTestPatient,'derivatives','ExploreASL','Population','qCBF_sub-Sub1_ASL_1.nii');

try
    % Load images
    imRefM0 = xASL_io_Nifti2Im(groundTruthM0File);
    imRefPerf = xASL_io_Nifti2Im(groundTruthPerfusionFile);
    imDerM0 = xASL_io_Nifti2Im(derivedM0File);
    imDerCBF = xASL_io_Nifti2Im(derivedCBFFile);
    
    % Test comparisons ... (not finished, work in progress)
    RMSE.M0 = xASL_ut_GetRMSE(imRefM0, imDerM0);
    RMSE.CBF = xASL_ut_GetRMSE(imRefPerf, imDerCBF);
    
    % Check RMSE values
    thresh.m0RMSE = 3;
    thresh.cbfRMSE = 3;
    if RMSE.M0>thresh.m0RMSE
        testCondition = false;
    end
    if RMSE.CBF>thresh.cbfRMSE
        testCondition = false;
    end
    
    % Mean values
    referenceMeanDerCBF = 8.8205;
    meanDerCBF = mean(imDerCBF(:));
    
    % Check that at least the mean value is somewhat consistent
    thresh.DiffMeanCBF = 0.1;
    absDiffMeanCBF = abs(referenceMeanDerCBF-meanDerCBF);
    if absDiffMeanCBF>thresh.DiffMeanCBF
        testCondition = false;
    end
    
catch ME
    warning(ME.identifier, '%s', ME.message);
    testCondition = false;
end


% Delete test data
xASL_delete(testPatientDestination,true)

% Clean-up
clearvars -except UnitTest TestRepository testCondition testTime testVersion

% Get test duration
UnitTest.tests(10).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(10).passed = testCondition;



%% Test run 11

% Give your individual subtest a name
UnitTest.tests(11).testname = 'DRO 2.3.0 (Pre-release version, multi-session BIDS to legacy)';

% Start the test
testTime = tic;

% Set-up DRO
subjectName = 'sub-001';
droTestPatientSource = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_3_0_pre_release');
droTestPatient = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_3_0_pre_release');
xASL_Copy(droTestPatientSource,fullfile(droTestPatient,'rawdata'),1);

% Set-up sessions
xASL_adm_CreateDir(fullfile(droTestPatient,'rawdata',subjectName,'ses-1'));
xASL_adm_CreateDir(fullfile(droTestPatient,'rawdata',subjectName,'ses-2'));
xASL_adm_CreateDir(fullfile(droTestPatient,'rawdata',subjectName,'ses-3'));
xASL_adm_CreateDir(fullfile(droTestPatient,'rawdata',subjectName,'ses-4'));
xASL_adm_CreateDir(fullfile(droTestPatient,'rawdata',subjectName,'ses-5'));

% Copy modalities

% Session one: only t1
xASL_Copy(fullfile(droTestPatient,'rawdata',subjectName,'anat'),fullfile(droTestPatient,'rawdata',subjectName,'ses-1','anat'));

% Session two: only asl
xASL_Copy(fullfile(droTestPatient,'rawdata',subjectName,'perf'),fullfile(droTestPatient,'rawdata',subjectName,'ses-2','perf'));

% Session three: t1 & asl
xASL_Copy(fullfile(droTestPatient,'rawdata',subjectName,'anat'),fullfile(droTestPatient,'rawdata',subjectName,'ses-3','anat'));
xASL_Copy(fullfile(droTestPatient,'rawdata',subjectName,'perf'),fullfile(droTestPatient,'rawdata',subjectName,'ses-3','perf'));

% Session four: t1, flair, & asl (using T2 for flair here)
xASL_Copy(fullfile(droTestPatient,'rawdata',subjectName,'anat'),fullfile(droTestPatient,'rawdata',subjectName,'ses-4','anat'));
xASL_Copy(fullfile(droTestPatient,'rawdata',subjectName,'perf'),fullfile(droTestPatient,'rawdata',subjectName,'ses-4','perf'));
xASL_Copy(fullfile(droTestPatient,'rawdata',subjectName,'ground_truth','sub-001_acq-002_T2map.nii.gz'),...
          fullfile(droTestPatient,'rawdata',subjectName,'ses-4','anat','sub-001_acq-002_FLAIR.nii.gz'));
xASL_Copy(fullfile(droTestPatient,'rawdata',subjectName,'ground_truth','sub-001_acq-002_T2map.json'),...
          fullfile(droTestPatient,'rawdata',subjectName,'ses-4','anat','sub-001_acq-002_FLAIR.json'));

% Session five: missing scans

% Delete templates
xASL_delete(fullfile(droTestPatient,'rawdata',subjectName,'anat'),true);
xASL_delete(fullfile(droTestPatient,'rawdata',subjectName,'perf'),true);
xASL_delete(fullfile(droTestPatient,'rawdata',subjectName,'ground_truth'),true);

% Fallback
testCondition = true;

% Test: BIDS2LEGACY
try
    [x] = ExploreASL_Master(droTestPatient,[0 0 0 1],0);
catch ME
    warning(ME.identifier, '%s', ME.message);
    testCondition = false;
    diary off;
end

% Test: Load data
try
    [x] = ExploreASL_Master(droTestPatient,0,0);
catch ME
    warning(ME.identifier, '%s', ME.message);
    testCondition = false;
    diary off;
end

% Test directories and files
testDirsAndFiles.session1 = fullfile(droTestPatient,'derivatives','ExploreASL','sub-001_1');
testDirsAndFiles.session2 = fullfile(droTestPatient,'derivatives','ExploreASL','sub-001_2');
testDirsAndFiles.session3 = fullfile(droTestPatient,'derivatives','ExploreASL','sub-001_3');
testDirsAndFiles.session4 = fullfile(droTestPatient,'derivatives','ExploreASL','sub-001_4');
testDirsAndFiles.session5 = fullfile(droTestPatient,'derivatives','ExploreASL','sub-001_5');
testDirsAndFiles.session1_t1w = fullfile(droTestPatient,'derivatives','ExploreASL','sub-001_1','T1.nii.gz');
testDirsAndFiles.session2_asl = fullfile(droTestPatient,'derivatives','ExploreASL','sub-001_2','ASL_1','ASL4D.nii.gz');
testDirsAndFiles.session2_m0 = fullfile(droTestPatient,'derivatives','ExploreASL','sub-001_2','ASL_1','M0.nii.gz');
testDirsAndFiles.session3_t1w = fullfile(droTestPatient,'derivatives','ExploreASL','sub-001_3','T1.nii.gz');
testDirsAndFiles.session3_asl = fullfile(droTestPatient,'derivatives','ExploreASL','sub-001_3','ASL_1','ASL4D.nii.gz');
testDirsAndFiles.session3_m0 = fullfile(droTestPatient,'derivatives','ExploreASL','sub-001_3','ASL_1','M0.nii.gz');
testDirsAndFiles.session4_t1w = fullfile(droTestPatient,'derivatives','ExploreASL','sub-001_4','T1.nii.gz');
testDirsAndFiles.session4_t1w = fullfile(droTestPatient,'derivatives','ExploreASL','sub-001_4','FLAIR.nii.gz');
testDirsAndFiles.session4_asl = fullfile(droTestPatient,'derivatives','ExploreASL','sub-001_4','ASL_1','ASL4D.nii.gz');
testDirsAndFiles.session4_m0 = fullfile(droTestPatient,'derivatives','ExploreASL','sub-001_4','ASL_1','M0.nii.gz');

% Define one or multiple test conditions here
if ~exist('x','var')
    testCondition = false;
else
    if ~isstruct(x)
        testCondition = false;
    end
end

% Iterate over test directories and files
fieldsTestDirsAndFiles = fieldnames(testDirsAndFiles);
for iField = 1:numel(fieldsTestDirsAndFiles)
    if ~xASL_exist(testDirsAndFiles.(fieldsTestDirsAndFiles{iField}),'file') && ~xASL_exist(testDirsAndFiles.(fieldsTestDirsAndFiles{iField}),'dir')
        testCondition = false;
    end
end

% Delete test data
xASL_delete(droTestPatient,true)

% Clean-up
clearvars -except UnitTest TestRepository testCondition testTime testVersion

% Get test duration
UnitTest.tests(11).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(11).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


%% Determine RMSE of two images
function RMSE = xASL_ut_GetRMSE(imageA, imageB)
    
    if isequal(size(imageA),size(imageB))
        % Calculate
        RMSE = sqrt(mean((imageA(:) - imageB(:)).^2))*2/sqrt(mean(abs(imageA(:)) + abs(imageB(:))).^2);
    else
        % Resample A to B
        [imageA] = xASL_im_ResampleIM(imageA, [1 0 0 0;0 1 0 0;0 0 1 0; 0 0 0 1], [1 0 0 0;0 1 0 0;0 0 1 0; 0 0 0 1], size(imageB));
        % Calculate
        RMSE = sqrt(mean((imageA(:) - imageB(:)).^2))*2/sqrt(mean(abs(imageA(:)) + abs(imageB(:))).^2);
    end
end

