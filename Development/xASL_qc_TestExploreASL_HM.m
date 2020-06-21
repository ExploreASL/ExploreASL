function [ResultsTable] = xASL_qc_TestExploreASL_HM(RunMethod, bTestSPM, bOverwrite) 
%xASL_qc_TestExploreASL_HM Run ExploreASL QC for test datasets

if nargin<1 || isempty(RunMethod)
    RunMethod = 1;
end
if nargin<2 || isempty(bTestSPM)
    bTestSPM = 1;
end
if nargin<3 || isempty(bOverwrite)
    bOverwrite = 1;
end

if ismac
    TestDirOrig = '/Users/henk/surfdrive/HolidayPics/ExploreASL_TestCases';
    TestDirDest = '/Users/henk/ExploreASL/ASL/TestCasesProcessed';
elseif ispc
    TestDirOrig = 'S:\gifmi\Projects\ExploreASL\ExploreASL_TestCases\ExploreASL_TestCases';
    TestDirDest = 'S:\gifmi\Projects\ExploreASL\ExploreASL_TestCases\ProcessedCases';
else
    % linux
     TestDirOrig = '/home/henk/ownCloud/HolidayPics/ExploreASL_TestCases';
     TestDirDest = '/home/henk/ExploreASL/ASL/ExploreASL_TestCasesProcessed';
end
   

[ResultsTable] = xASL_qc_TestExploreASL(TestDirOrig, TestDirDest, RunMethod, bTestSPM, [], [], [], bOverwrite)    


end

