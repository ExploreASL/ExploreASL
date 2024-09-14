function UnitTest = xASL_ut_CheckSubtests(UnitTest)
%xASL_ut_CheckSubtests Check individual subtests
%
% INPUT:        UnitTest  - Test structure
%
% OUTPUT:       UnitTest  - Test structure
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Check if an individual subtest failed.
%
% EXAMPLE:      UnitTest = xASL_ut_CheckSubtests(UnitTest);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

    % Check if an individual subtest failed
    UnitTest.passed = true;
    for it = 1:numel(UnitTest.tests)
        if ~UnitTest.tests(it).passed
            UnitTest.passed = false;
        end
    end
    
end
