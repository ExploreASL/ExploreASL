function UnitTest = xASL_ut_function_xASL_bids_VendorFieldCheck(TestRepository)
%xASL_ut_function_xASL_bids_VendorFieldCheck Individual unit test for xASL_bids_VendorFieldCheck
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_bids_VendorFieldCheck(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Test GE JSON';

% Start the test
testTime = tic;

% Run your test here
jsonIn.Manufacturer = 'GE';
jsonIn.CoilString = 'ThisIsACoilString';
jsonIn.EffectiveEchoSpacing = -1;
jsonIn.TotalReadoutTime = -1;
jsonIn.NumberSegments = 1;
jsonIn.PhaseEncodingAxis = 'i';
jsonOut = xASL_bids_VendorFieldCheck(jsonIn);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~isstruct(jsonOut)
    testCondition = false;
end
% Check if all fields exist
if ~isfield(jsonOut,'ReceiveCoilName') || ...
   ~isfield(jsonOut,'NumberShots') || ...
   ~isfield(jsonOut,'PhaseEncodingDirection') || ...
   ~isfield(jsonOut,'Manufacturer') || ...
   ~isfield(jsonOut,'EffectiveEchoSpacing') || ...
   ~isfield(jsonOut,'TotalReadoutTime')
    testCondition = false;
end
% Check field content
if isfield(jsonOut,'ReceiveCoilName') && ...
   isfield(jsonOut,'NumberShots') && ...
   isfield(jsonOut,'PhaseEncodingDirection') && ...
   isfield(jsonOut,'Manufacturer') && ...
   isfield(jsonOut,'EffectiveEchoSpacing') && ...
   isfield(jsonOut,'TotalReadoutTime')
    % Check each field
    if ~strcmp(jsonOut.ReceiveCoilName,'ThisIsACoilString')
        testCondition = false;
    end
    if ~(jsonOut.NumberShots==1)
        testCondition = false;
    end
    if ~strcmp(jsonOut.PhaseEncodingDirection,'i')
        testCondition = false;
    end
    if ~strcmp(jsonOut.Manufacturer,'GE')
        testCondition = false;
    end
    if ~(jsonOut.EffectiveEchoSpacing==1)
        testCondition = false;
    end
    if ~(jsonOut.TotalReadoutTime==1)
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


