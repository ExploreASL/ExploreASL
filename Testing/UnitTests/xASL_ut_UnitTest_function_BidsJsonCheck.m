function UnitTest = xASL_ut_UnitTest_function_BidsJsonCheck(TestRepository)
%xASL_ut_UnitTest_function_BidsJsonCheck Individual unit test for xASL_bids_JsonCheck
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
% EXAMPLE:      UnitTests(1) = xASL_ut_UnitTest_function_BidsJsonCheck(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL

%% Initialize test structure

% Insert test name here
UnitTest.name = 'xASL_bids_JsonCheck';

% Define whether you are testing a module, submodule or function
UnitTest.unit = 'Function';

%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Check empty ASL JSON';

% Start the test
testTime = tic;

% Run your test here
jsonIn = struct;
jsonOut = xASL_bids_JsonCheck(jsonIn,'ASL');

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~isempty(jsonOut)
    testCondition = false;
end

% Clean up
clear jsonIn jsonOut

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% Test run 2

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Check non-empty ASL JSON';

% Start the test
testTime = tic;

% Run your test here
jsonIn = struct;
jsonIn.TotalAcquiredPairs = 1;
jsonIn.RepetitionTimePreparation = 1;
jsonIn.EchoTime = 1;
jsonIn.MagneticFieldStrength = 1;
jsonIn.MRAcquisitionType = '3D';
jsonIn.ArterialSpinLabelingType = 'PCASL';
jsonIn.PostLabelingDelay = 1;
jsonIn.BackgroundSuppression = 0;
jsonIn.M0Type = 'Separate';
jsonIn.LabelingDuration = 1;
jsonOut = xASL_bids_JsonCheck(jsonIn,'ASL');

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~isstruct(jsonOut)
    testCondition = false;
end
% Check fields
if ~isfield(jsonOut,'MagneticFieldStrength') || ...
   ~isfield(jsonOut,'MRAcquisitionType') || ...
   ~isfield(jsonOut,'EchoTime') || ...
   ~isfield(jsonOut,'RepetitionTimePreparation') || ...
   ~isfield(jsonOut,'ArterialSpinLabelingType') || ...
   ~isfield(jsonOut,'PostLabelingDelay') || ...
   ~isfield(jsonOut,'BackgroundSuppression') || ...
   ~isfield(jsonOut,'M0Type') || ...
   ~isfield(jsonOut,'TotalAcquiredPairs') || ...
   ~isfield(jsonOut,'LabelingDuration')
    testCondition = false;
end

% Clean up
clear jsonIn jsonOut

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;


%% Test run 3

% Give your individual subtest a name
UnitTest.tests(3).testname = 'Check non-empty M0 JSON';

% Start the test
testTime = tic;

% Run your test here
jsonIn = struct;
jsonIn.RepetitionTimePreparation = 1;
jsonIn.EchoTime = 1;
jsonIn.MagneticFieldStrength = 1;
jsonIn.MRAcquisitionType = '3D';
jsonIn.ArterialSpinLabelingType = 'PCASL';
jsonIn.PostLabelingDelay = 1;
jsonIn.BackgroundSuppression = 0;
jsonIn.M0Type = 'Separate';
jsonIn.LabelingDuration = 1;
jsonOut = xASL_bids_JsonCheck(jsonIn,'M0');

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~isstruct(jsonOut)
    testCondition = false;
end
% Check fields
if ~isfield(jsonOut,'MagneticFieldStrength') || ...
   ~isfield(jsonOut,'MRAcquisitionType') || ...
   ~isfield(jsonOut,'EchoTime') || ...
   ~isfield(jsonOut,'RepetitionTimePreparation') || ...
   ~isfield(jsonOut,'ArterialSpinLabelingType') || ...
   ~isfield(jsonOut,'PostLabelingDelay') || ...
   ~isfield(jsonOut,'BackgroundSuppression') || ...
   ~isfield(jsonOut,'M0Type') || ...
   ~isfield(jsonOut,'LabelingDuration')
    testCondition = false;
end

% Clean up
clear jsonIn jsonOut

% Get test duration
UnitTest.tests(3).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(3).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


