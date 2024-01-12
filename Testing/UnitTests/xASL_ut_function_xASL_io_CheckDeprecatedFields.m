function UnitTest = xASL_ut_function_xASL_io_CheckDeprecatedFields(TestRepository)
%xASL_ut_function_xASL_io_CheckDeprecatedFields Individual unit test for xASL_io_CheckDeprecatedFields
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_io_CheckDeprecatedFields(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2023 ExploreASL


%% Test run 1

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Deprecated field renaming';

% Start the test
testTime = tic;

% Prepare the x-struct with deprecated names
x.RepetitionTime = 4400;
x.EchoTime = 14.3660;
x.NumberOfTemporalPositions = 140;
x.RescaleSlope = 10.8730;
x.RescaleSlopeOriginal = 10.8730;
x.MRScaleSlope = 0.0095;
x.RescaleIntercept = 0;
x.AcquisitionTime = 1.3365e+05;
x.SavePWI4D = 1;
x.Segment_SPM12 = 1;
x.M0_conventionalProcessing = 1;
x.M0 = 'separate_scan';
x.Sequence = '2D_EPI';
x.Lambda = 0.9000;
x.T2 = 165;
x.TissueT1 = 1300;
x.nCompartments = 2;
x.readoutDim = '2D';
x.Vendor = 'Philips';
x.Q.LabelingType = 'CASL';
x.Q.LabelingDuration = 1800;
x.Q.BloodT1 = 1818;
x.Q.BackgroundSuppressionNumberPulses = 2;
x.Q.Initial_PLD = 1800;

% Run the conversion
xResult = xASL_io_CheckDeprecatedFieldsX(x);

% Create the reference file
xReference.RepetitionTime = 4400;
xReference.Q.EchoTime = 14.3660;
xReference.NumberOfTemporalPositions = 140;
xReference.RescaleSlope = 10.8730;
xReference.RescaleSlopeOriginal = 10.8730;
xReference.MRScaleSlope = 0.0095;
xReference.RescaleIntercept = 0;
xReference.AcquisitionTime = 133650;
xReference.SavePWI4D = 1;
xReference.T2 = 165;
xReference.Q.LabelingType = 'CASL';
xReference.Q.LabelingDuration = 1800;
xReference.Q.BloodT1 = 1818;
xReference.Q.BackgroundSuppressionNumberPulses = 2;
xReference.Q.Initial_PLD = 1800;
xReference.Q.Sequence = '2D_EPI';
xReference.Q.Vendor = 'Philips';
xReference.Q.readoutDim = '2D';
xReference.Q.Lambda = 0.9000;
xReference.Q.TissueT1 = 1300;
xReference.Q.nCompartments = 2;
xReference.Q.M0 = 'separate_scan';
xReference.modules.structural.bSegmentSPM12 = 1;
xReference.modules.asl.M0_conventionalProcessing = 1;

% Compare the result with the reference
testCondition = true; % Fallback

if isempty(xResult)
	testCondition = false;
end

if ~isequal(xResult, xReference)
	testCondition = false;
end

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;

%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end
