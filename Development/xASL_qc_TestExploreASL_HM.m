function [ResultsTable] = xASL_qc_TestExploreASL_HM
%xASL_qc_TestExploreASL_HM Run ExploreASL QC for test datasets


if ismac
    TestDirOrig='/Users/henk/surfdrive/HolidayPics/ExploreASL_TestCases';
    TestDirDest='/Users/henk/ExploreASL/ASL/TestCasesProcessed';
elseif ispc
    TestDirOrig='/Users/henk/surfdrive/HolidayPics/ExploreASL_TestCases';
    TestDirDest='/Users/henk/ExploreASL/ASL/TestCasesProcessed';
else
    % linux
    [StatusIs, UserName] = system('whoami','-echo');

     TestDirOrig='/home/henk/ownCloud/HolidayPics/ExploreASL_TestCases';
     TestDirDest='/home/henk/ExploreASL/ASL/ExploreASL_TestCasesProcessed';
end
    
    
    
[ResultsTable] = xASL_qc_TestExploreASL(TestDirOrig, TestDirDest, 1, 1);


end

