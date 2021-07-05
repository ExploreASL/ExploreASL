function UnitTest = xASL_ut_function_xASL_adm_OrderFields(TestRepository)
%xASL_ut_function_xASL_adm_OrderFields Individual unit test for xASL_adm_OrderFields
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_adm_OrderFields(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Extract integer example';

% Start the test
testTime = tic;

% Run your test here
inStruct.A = 'test';
inStruct.D = 'test';
inStruct.C = 'test';
inStruct.B = 'test';
orderStruct.A = 'test';
orderStruct.B = 'test';
orderStruct.C = 'test';
orderStruct.D = 'test';
outStruct = xASL_adm_OrderFields(inStruct,orderStruct);
fieldNamesOutStruct = fieldnames(outStruct);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~strcmp(fieldNamesOutStruct(1),'A') || ...
   ~strcmp(fieldNamesOutStruct(2),'B') || ...
   ~strcmp(fieldNamesOutStruct(3),'C') || ...
   ~strcmp(fieldNamesOutStruct(4),'D')
    testCondition = false;
end

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


