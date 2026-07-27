function UnitTest = xASL_ut_function_xASL_adm_MatlabVersionYear(TestRepository)
%xASL_ut_function_xASL_adm_MatlabVersionYear Unit test for xASL_adm_MatlabVersionYear
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_adm_MatlabVersionYear(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________


%% ---- Group A: Return type and validity ---------------------------------

% Subtest 1: Return value is a numeric scalar
UnitTest.tests(1).testname = 'Return value is a numeric scalar';
testTime = tic;
matlabVersionYear = xASL_adm_MatlabVersionYear;
UnitTest.tests(1).passed = isnumeric(matlabVersionYear) && isscalar(matlabVersionYear);
UnitTest.tests(1).duration = toc(testTime);

% Subtest 2: MATLAB release year is plausible (>=2019)
UnitTest.tests(2).testname = 'MATLAB release year is >= 2019';
testTime = tic;
matlabVersionYear = xASL_adm_MatlabVersionYear;
UnitTest.tests(2).passed = matlabVersionYear >= 2019;
UnitTest.tests(2).duration = toc(testTime);

% Subtest 3: MATLAB release year is not unrealistically far in the future
UnitTest.tests(3).testname = 'MATLAB release year is plausible';
testTime = tic;
matlabVersionYear = xASL_adm_MatlabVersionYear;
currentYear = year(datetime('today'));
UnitTest.tests(3).passed = matlabVersionYear <= (currentYear + 1);
UnitTest.tests(3).duration = toc(testTime);


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end