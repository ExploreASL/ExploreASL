function UnitTest = xASL_ut_function_xASL_bids_CheckDatasetDescription(TestRepository)
%xASL_ut_function_xASL_bids_CheckDatasetDescription Individual unit test for xASL_bids_CheckDatasetDescription
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_bids_CheckDatasetDescription(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Correct version, different version, non-xASL import, missing field';

% Start the test
testTime = tic;

% Run your test here
datasetDescription.Name = 'TestPatient';
datasetDescription.BIDSVersion = '1.6.0';
datasetDescription.DatasetType = 'raw';
datasetDescription.HowToAcknowledge = 'Please cite this paper: https://www.ncbi.nlm.nih.gov/pubmed/001012092119281';

% Example 1: correct version
datasetDescription.Acknowledgements = 'Imported with ExploreASL 1.8.0_BETA';
xVersion = '1.8.0_BETA';
[bImportedExploreASL1, bImportedSameVersion1, versionExploreASLBIDS1, bImportedBETA1] = xASL_bids_CheckDatasetDescription(datasetDescription,xVersion);

% Example 2: different version
datasetDescription.Acknowledgements = 'Imported with ExploreASL 1.7.0';
xVersion = '1.5.0_BETA';
[bImportedExploreASL2, bImportedSameVersion2, versionExploreASLBIDS2, bImportedBETA2] = xASL_bids_CheckDatasetDescription(datasetDescription,xVersion);

% Example 3: non-xASL import
datasetDescription.Acknowledgements = 'Not with ExploreASL';
xVersion = '1.9.0';
[bImportedExploreASL3, bImportedSameVersion3, versionExploreASLBIDS3, bImportedBETA3] = xASL_bids_CheckDatasetDescription(datasetDescription,xVersion);

% Example 4: missing field
datasetDescription = rmfield(datasetDescription,'Acknowledgements');
xVersion = '1.3.0';
[bImportedExploreASL4, bImportedSameVersion4, versionExploreASLBIDS4, bImportedBETA4] = xASL_bids_CheckDatasetDescription(datasetDescription,xVersion);


% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here

% Check file types first
if ~islogical(bImportedExploreASL1) || ~islogical(bImportedExploreASL2) || ~islogical(bImportedExploreASL3) || ~islogical(bImportedExploreASL4) 
    testCondition = false;
end
if ~islogical(bImportedSameVersion1) || ~islogical(bImportedSameVersion2) || ~islogical(bImportedSameVersion3) || ~islogical(bImportedSameVersion4) 
    testCondition = false;
end
if ~ischar(versionExploreASLBIDS1) || ~ischar(versionExploreASLBIDS2) || ~isempty(versionExploreASLBIDS3) || ~isempty(versionExploreASLBIDS4) 
    testCondition = false;
end
if ~islogical(bImportedBETA1) || ~islogical(bImportedBETA2) || ~islogical(bImportedBETA3) || ~islogical(bImportedBETA4) 
    testCondition = false;
end

% Check content
if ~strcmp(versionExploreASLBIDS1,'1.8.0') || ~strcmp(versionExploreASLBIDS2,'1.7.0')
    testCondition = false;
end
if ~bImportedBETA1 || bImportedBETA2 || bImportedBETA3 || bImportedBETA4
    testCondition = false;
end
if ~bImportedExploreASL1 || ~bImportedExploreASL2 || bImportedExploreASL3 || bImportedExploreASL4
    testCondition = false;
end
if ~bImportedSameVersion1 || bImportedSameVersion2 || bImportedSameVersion3 || bImportedSameVersion4
    testCondition = false;
end

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


