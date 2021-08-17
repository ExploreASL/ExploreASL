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
datasetDescription.Acknowledgements = 'Imported with ExploreASL 1.8.0_BETA.';
xVersion = '1.8.0_BETA';
[thisTest.bImportedExploreASL1, thisTest.bImportedSameVersion1, thisTest.versionExploreASLBIDS1, thisTest.bImportedBETA1] = xASL_bids_CheckDatasetDescription(datasetDescription,xVersion);

% Example 2: different version
datasetDescription.Acknowledgements = 'Imported with ExploreASL 1.7.0.';
xVersion = '1.5.0_BETA';
[thisTest.bImportedExploreASL2, thisTest.bImportedSameVersion2, thisTest.versionExploreASLBIDS2, thisTest.bImportedBETA2] = xASL_bids_CheckDatasetDescription(datasetDescription,xVersion);

% Example 3: non-xASL import
datasetDescription.Acknowledgements = 'Not with ExploreASL.';
xVersion = '1.9.0';
[thisTest.bImportedExploreASL3, thisTest.bImportedSameVersion3, thisTest.versionExploreASLBIDS3, thisTest.bImportedBETA3] = xASL_bids_CheckDatasetDescription(datasetDescription,xVersion);

% Example 4: missing field
datasetDescription = rmfield(datasetDescription,'Acknowledgements');
xVersion = '1.3.0';
[thisTest.bImportedExploreASL4, thisTest.bImportedSameVersion4, thisTest.versionExploreASLBIDS4, thisTest.bImportedBETA4] = xASL_bids_CheckDatasetDescription(datasetDescription,xVersion);


% Define one or multiple test conditions here
testCondition = true;

% Define one or multiple test conditions here

% Check file types first
if ~islogical(thisTest.bImportedExploreASL1) ...
        || ~islogical(thisTest.bImportedExploreASL2) ...
        || ~islogical(thisTest.bImportedExploreASL3) ...
        || ~islogical(thisTest.bImportedExploreASL4) 
    testCondition = false;
end
if ~islogical(thisTest.bImportedSameVersion1) ...
        || ~islogical(thisTest.bImportedSameVersion2) ...
        || ~islogical(thisTest.bImportedSameVersion3) ...
        || ~islogical(thisTest.bImportedSameVersion4) 
    testCondition = false;
end
if ~ischar(thisTest.versionExploreASLBIDS1) ...
        || ~ischar(thisTest.versionExploreASLBIDS2) ...
        || ~isempty(thisTest.versionExploreASLBIDS3) ...
        || ~isempty(thisTest.versionExploreASLBIDS4) 
    testCondition = false;
end
if ~islogical(thisTest.bImportedBETA1) ...
        || ~islogical(thisTest.bImportedBETA2) ...
        || ~islogical(thisTest.bImportedBETA3) ...
        || ~islogical(thisTest.bImportedBETA4) 
    testCondition = false;
end

% Check content
if ~strcmp(thisTest.versionExploreASLBIDS1,'1.8.0_BETA') || ~strcmp(thisTest.versionExploreASLBIDS2,'1.7.0')
    testCondition = false;
end
if ~thisTest.bImportedBETA1 || thisTest.bImportedBETA2 || thisTest.bImportedBETA3 || thisTest.bImportedBETA4
    testCondition = false;
end
if ~thisTest.bImportedExploreASL1 || ~thisTest.bImportedExploreASL2 || thisTest.bImportedExploreASL3 || thisTest.bImportedExploreASL4
    testCondition = false;
end
if ~thisTest.bImportedSameVersion1 || thisTest.bImportedSameVersion2 || thisTest.bImportedSameVersion3 || thisTest.bImportedSameVersion4
    testCondition = false;
end

% Clean-up
clear datasetDescription xVersion thisTest

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;


%% Test run 2

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Random text before and after ExploreASL statement';

% Start the test
testTime = tic;

% Run your test here
datasetDescription.Name = 'TestPatient';
datasetDescription.BIDSVersion = '1.6.0';
datasetDescription.DatasetType = 'raw';
datasetDescription.HowToAcknowledge = 'Please cite this paper: https://www.ncbi.nlm.nih.gov/pubmed/001012092119281';

% Example 1: random text before and after with beta
datasetDescription.Acknowledgements = 'ThisIsSomeRandomText.Imported with ExploreASL 1.8.0_BETA.ThisIsEvenMoreText.';
xVersion = '1.8.0_BETA';
[thisTest.bImportedExploreASL1, thisTest.bImportedSameVersion1, thisTest.versionExploreASLBIDS1, thisTest.bImportedBETA1] = xASL_bids_CheckDatasetDescription(datasetDescription,xVersion);

% Example 2: random text before and after with multi-digit version number and without beta
datasetDescription.Acknowledgements = 'ThisIsSomeRandomText.Imported with ExploreASL 1.10.0.ThisIsEvenMoreText.';
xVersion = '1.8.0_BETA';
[thisTest.bImportedExploreASL2, thisTest.bImportedSameVersion2, thisTest.versionExploreASLBIDS2, thisTest.bImportedBETA2] = xASL_bids_CheckDatasetDescription(datasetDescription,xVersion);

% Define one or multiple test conditions here
testCondition = true;

% Check file types first
if ~islogical(thisTest.bImportedExploreASL1) ...
        || ~islogical(thisTest.bImportedExploreASL2) 
    testCondition = false;
end
if ~strcmp(thisTest.versionExploreASLBIDS1,'1.8.0_BETA') ...
        || ~strcmp(thisTest.versionExploreASLBIDS2,'1.10.0')
    testCondition = false;
end
if ~thisTest.bImportedBETA1 || thisTest.bImportedBETA2
    testCondition = false;
end
if ~thisTest.bImportedSameVersion1 || thisTest.bImportedSameVersion2
    testCondition = false;
end

% Clean-up
clear datasetDescription xVersion thisTest

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


