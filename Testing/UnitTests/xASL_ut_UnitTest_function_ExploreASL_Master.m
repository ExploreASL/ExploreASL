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
if isfield(x, 'opts')
    if ~isfield(x.opts,'DataParPath'),     testCondition = false;      end
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
    if isfield(x.opts,'DataParPath')
        if ~isempty(x.opts.DataParPath) || ~ischar(x.opts.DataParPath)
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
if isfield(x, 'opts')
    if ~isfield(x.opts,'DataParPath'),     testCondition = false;      end
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
    if isfield(x.opts,'DataParPath')
        if ~isempty(x.opts.DataParPath) || ~ischar(x.opts.DataParPath)
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
if isfield(x, 'opts')
    if ~isfield(x.opts,'DataParPath'),     testCondition = false;      end
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
    if isfield(x.opts,'DataParPath')
        if ~isempty(x.opts.DataParPath) || ~ischar(x.opts.DataParPath)
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
    warning('%s', ME.message);
    testCondition = false;
    diary off;
end

% Define one or multiple test conditions here
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

% Give your individual subtest a name
UnitTest.tests(5).testname = 'DRO 2.2.0 (NII2BIDS on rawdata)';

% Start the test
testTime = tic;

% Set-up DRO
droTestPatientSource = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_2_0');
droTestPatient = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_2_0');
droSubject = 'sub-Sub1'; % DRO subject
xASL_Copy(droTestPatientSource,fullfile(droTestPatient,'rawdata',droSubject),1);
xASL_bids_DRO2BIDS(droTestPatient); % Prepare DRO

% Fallback
testCondition = true;

% Read test files
try
    [x] = ExploreASL_Master(testPatientDestination,[0 1 0 0],0,0,1,1);
catch ME
    warning('%s', ME.message);
    testCondition = false;
    diary off;
end

% Define one or multiple test conditions here

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

% Get test duration
UnitTest.tests(5).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(5).passed = testCondition;


%% Test run 6

% Give your individual subtest a name
UnitTest.tests(6).testname = 'DRO 2.2.0 (ANONYMIZE, BIDS2Legacy)';

% Start the test
testTime = tic;

% Set-up DRO
droTestPatientSource = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_2_0');
droTestPatient = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_2_0');
droSubject = 'sub-Sub1'; % DRO subject
xASL_Copy(droTestPatientSource,fullfile(droTestPatient,'rawdata',droSubject),1);
xASL_bids_DRO2BIDS(droTestPatient); % Prepare DRO

% Fallback
testCondition = true;

% Read test files
try
    [x] = ExploreASL_Master(testPatientDestination,[0 1 1 1],0,0,1,1);
catch ME
    warning('%s', ME.message);
    testCondition = false;
    diary off;
end

% Define one or multiple test conditions here

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

% Get test duration
UnitTest.tests(6).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(6).passed = testCondition;


%% Test run 7

% Give your individual subtest a name
UnitTest.tests(7).testname = 'DRO 2.2.0 (ANONYMIZE, BIDS2Legacy with specific dataPar.json)';

% Start the test
testTime = tic;

% Set-up DRO
droTestPatientSource = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_2_0');
droTestPatient = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_2_0');
droSubject = 'sub-Sub1'; % DRO subject
xASL_Copy(droTestPatientSource,fullfile(droTestPatient,'rawdata',droSubject),1);
xASL_bids_DRO2BIDS(droTestPatient); % Prepare DRO
% Create dataPar.json
dataParStruct.x.Quality = 0;
dataParStruct.x.S.Atlases = {'TotalGM','DeepWM','Hammers','HOcort_CONN','HOsub_CONN','Mindboggle_OASIS_DKT31_CMA'};
spm_jsonwrite(fullfile(droTestPatient,'dataPar.json'),dataParStruct);

% Fallback
testCondition = true;

% Read test files
try
    [x] = ExploreASL_Master(testPatientDestination,[0 0 1 1],0,0,1,1);
catch ME
    warning('%s', ME.message);
    testCondition = false;
    diary off;
end

% Define one or multiple test conditions here

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

% Get test duration
UnitTest.tests(7).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(7).passed = testCondition;


%% Test run 8

% Give your individual subtest a name
UnitTest.tests(8).testname = 'DRO 2.2.0 (full pipeline, NII->Results)';

% Start the test
testTime = tic;

% Set-up DRO
droTestPatientSource = fullfile(TestRepository,'UnitTesting','dro_files','test_patient_2_2_0');
droTestPatient = fullfile(TestRepository,'UnitTesting','working_directory','test_patient_2_2_0');
droSubject = 'sub-Sub1'; % DRO subject
xASL_Copy(droTestPatientSource,fullfile(droTestPatient,'rawdata',droSubject),1);
xASL_bids_DRO2BIDS(droTestPatient,droSubject,false); % Prepare DRO, keep ground truth data for comparison

% Fallback
testCondition = true;

% Read test files
try
    [x] = ExploreASL_Master(testPatientDestination,[0 1 1 1],1,0,1,1);
catch ME
    warning('%s', ME.message);
    testCondition = false;
    diary off;
end

% Define one or multiple test conditions here

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
    
    % ...
    
catch ME
    warning('%s', ME.message);
    testCondition = false;
end


% Delete test data
xASL_delete(testPatientDestination,true)

% Get test duration
UnitTest.tests(8).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(8).passed = testCondition;


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

