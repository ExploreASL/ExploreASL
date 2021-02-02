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

    % Check if an individual subtest failed
    UnitTest.passed = true;
    for it = 1:numel(UnitTest.tests)
        if ~UnitTest.tests(it).passed
            UnitTest.passed = false;
        end
    end
    
end

