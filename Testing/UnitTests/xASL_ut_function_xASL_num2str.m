function UnitTest = xASL_ut_function_xASL_num2str(TestRepository)
%xASL_ut_function_xASL_num2str Individual unit test for xASL_num2str
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_str2num(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Row vector, column vector, scalar number';

% Start the test
testTime = tic;

% Run your test here
outA = xASL_num2str([1, 2, 3]');
outB = xASL_num2str([1, 2, 3]);
outC = xASL_num2str(1);

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~ischar(outA) || ~ischar(outB) || ~ischar(outC)
    testCondition = false;
end
if isempty(regexp(outA,'1,2,3', 'once'))
    testCondition = false;
end
if isempty(regexp(outB,'1,2,3', 'once'))
    testCondition = false;
end
if ~strcmp(outC,'1')
    testCondition = false;
end

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% Test run 2

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Examples from header';

% Start the test
testTime = tic;

% Run your test here

% DataOut1 = '10.5798'
DataOut1 = xASL_num2str(10.5798);

% DataOut2 = '11'
DataOut2 = xASL_num2str(10.5798, 2);

% DataOut3 = '1,15'
DataOut3 = xASL_num2str([1;15], 2,1);

% DataOut4 = '1,,15,,2,,13'
DataOut4 = xASL_num2str([1,15;2,13], 2,1,',,');

% Get exactly 5 digits after the comma: DataOut5 = '123.45679'
DataOut5 = xASL_num2str(123.456789, '%.5f');

% Automatic mode (remove trailing zeros): DataOut6 = '1.23456789'
DataOut6 = xASL_num2str(1.23456789000, 'auto');

% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here
if ~ischar(DataOut1) || ~ischar(DataOut2) || ~ischar(DataOut3) || ...
   ~ischar(DataOut4) || ~ischar(DataOut5) || ~ischar(DataOut6)   
    testCondition = false;
end
if ~strcmp(DataOut1,'10.5798')
    testCondition = false;
end
if ~strcmp(DataOut2,'11')
    testCondition = false;
end
if ~strcmp(DataOut3,'1,15')
    testCondition = false;
end
if ~strcmp(DataOut4,'1,,15,,2,,13')
    testCondition = false;
end
if ~strcmp(DataOut5,'123.45679')
    testCondition = false;
end
if ~strcmp(DataOut6,'1.23456789')
    testCondition = false;
end

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


