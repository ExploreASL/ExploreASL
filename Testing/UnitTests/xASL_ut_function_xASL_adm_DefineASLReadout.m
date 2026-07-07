function UnitTest = xASL_ut_function_xASL_adm_DefineASLReadout(TestRepository)
%xASL_ut_function_xASL_adm_DefineASLReadout Individual unit test for xASL_adm_DefineASLReadout
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
% DESCRIPTION:  Should be run using xASL_ut_UnitTesting. Covers the tiered
%               PulseSequenceType classification (vendor-specific names +
%               Damerau-Levenshtein fuzzy fallback), the Manufacturer->Vendor
%               inference with fuzzy fallback, and the PulseSequenceType vs
%               MRAcquisitionType consistency validation.
%
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_adm_DefineASLReadout(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________



%% Suppress warnings during tests (warnings are tested implicitly by
% verifying the corrected field values; we don't want warning text to
% pollute the test output or to be mistaken for errors).
warnState = warning('off', 'all');



%% Test run 1: Motivating bug — typo "3D_SPRIAL" + Manufacturer=Siemens -> spiral, Vendor=Siemens

UnitTest.tests(1).testname = 'Typo 3D_SPRIAL + Manufacturer Siemens -> spiral, Vendor Siemens';

testTime = tic;

xQ = struct;
xQ.PulseSequenceType = '3D_SPRIAL';
xQ.Manufacturer      = 'Siemens';
xQ.MRAcquisitionType = '3D';
xQ = xASL_adm_DefineASLReadout(xQ, false);

testCondition = true;
if ~isfield(xQ, 'PulseSequenceType') || ~strcmpi(xQ.PulseSequenceType, 'spiral')
    testCondition = false;
end
if ~isfield(xQ, 'Vendor') || ~strcmpi(xQ.Vendor, 'Siemens')
    testCondition = false;
end

UnitTest.tests(1).duration = toc(testTime);
UnitTest.tests(1).passed = testCondition;


%% Test run 2: Siemens vendor-specific "tgse" -> GRASE

UnitTest.tests(2).testname = 'Siemens tgse -> GRASE';

testTime = tic;

xQ = struct;
xQ.PulseSequenceType = 'tgse';
xQ.Vendor            = 'Siemens';
xQ.MRAcquisitionType = '3D';
xQ = xASL_adm_DefineASLReadout(xQ, false);

testCondition = true;
if ~isfield(xQ, 'PulseSequenceType') || ~strcmpi(xQ.PulseSequenceType, 'GRASE')
    testCondition = false;
end

UnitTest.tests(2).duration = toc(testTime);
UnitTest.tests(2).passed = testCondition;


%% Test run 3: Siemens EPI family "ep2d_asl" -> EPI

UnitTest.tests(3).testname = 'Siemens ep2d_asl -> EPI';

testTime = tic;

xQ = struct;
xQ.PulseSequenceType = 'ep2d_asl';
xQ.Vendor            = 'Siemens';
xQ.MRAcquisitionType = '2D';
xQ = xASL_adm_DefineASLReadout(xQ, false);

testCondition = true;
if ~isfield(xQ, 'PulseSequenceType') || ~strcmpi(xQ.PulseSequenceType, 'EPI')
    testCondition = false;
end

UnitTest.tests(3).duration = toc(testTime);
UnitTest.tests(3).passed = testCondition;


%% Test run 4: Philips "FFE-EPI" -> EPI

UnitTest.tests(4).testname = 'Philips FFE-EPI -> EPI';

testTime = tic;

xQ = struct;
xQ.PulseSequenceType = 'FFE-EPI';
xQ.Vendor            = 'Philips';
xQ.MRAcquisitionType = '2D';
xQ = xASL_adm_DefineASLReadout(xQ, false);

testCondition = true;
if ~isfield(xQ, 'PulseSequenceType') || ~strcmpi(xQ.PulseSequenceType, 'EPI')
    testCondition = false;
end

UnitTest.tests(4).duration = toc(testTime);
UnitTest.tests(4).passed = testCondition;


%% Test run 5: Manufacturer present should protect Vendor from spiral->GE heuristic

UnitTest.tests(5).testname = 'spiral + Manufacturer Siemens preserves Vendor=Siemens (not GE)';

testTime = tic;

xQ = struct;
xQ.PulseSequenceType = 'spiral';
xQ.Manufacturer      = 'Siemens';
xQ.MRAcquisitionType = '3D';
xQ = xASL_adm_DefineASLReadout(xQ, false);

testCondition = true;
if ~isfield(xQ, 'Vendor') || ~strcmpi(xQ.Vendor, 'Siemens')
    testCondition = false;
end
if strcmpi(xQ.Vendor, 'GE')
    testCondition = false; % Vendor must NOT be auto-set to GE when Manufacturer is present
end

UnitTest.tests(5).duration = toc(testTime);
UnitTest.tests(5).passed = testCondition;


%% Test run 6: Manufacturer typo "Seimens" fuzzy-corrects to Vendor Siemens

UnitTest.tests(6).testname = 'Manufacturer typo Seimens -> Vendor Siemens via fuzzy match';

testTime = tic;

xQ = struct;
xQ.Manufacturer      = 'Seimens';
xQ.MRAcquisitionType = '3D';
xQ = xASL_adm_DefineASLReadout(xQ, false);

testCondition = true;
if ~isfield(xQ, 'Vendor') || ~strcmpi(xQ.Vendor, 'Siemens')
    testCondition = false;
end

UnitTest.tests(6).duration = toc(testTime);
UnitTest.tests(6).passed = testCondition;


%% Test run 7: spiral + neither Vendor nor Manufacturer -> fallback GE heuristic

UnitTest.tests(7).testname = 'spiral alone (no Vendor/Manufacturer) -> fallback Vendor GE';

testTime = tic;

xQ = struct;
xQ.PulseSequenceType = 'spiral';
xQ.MRAcquisitionType = '3D';
xQ = xASL_adm_DefineASLReadout(xQ, false);

testCondition = true;
if ~isfield(xQ, 'Vendor') || ~strcmpi(xQ.Vendor, 'GE')
    testCondition = false;
end

UnitTest.tests(7).duration = toc(testTime);
UnitTest.tests(7).passed = testCondition;


%% Test run 8: Consistency check — 3D + EPI on non-DRO vendor warns (no field change, just warn)

UnitTest.tests(8).testname = '3D + EPI (non-DRO) accepted; value preserved as EPI';

testTime = tic;

xQ = struct;
xQ.PulseSequenceType = '3D_EPI';
xQ.Vendor            = 'Siemens';
xQ.MRAcquisitionType = '3D';
xQ = xASL_adm_DefineASLReadout(xQ, false);

testCondition = true;
% The classification should map "3D_EPI" to "EPI"; the consistency check
% only warns, it must not erase the value.
if ~isfield(xQ, 'PulseSequenceType') || ~strcmpi(xQ.PulseSequenceType, 'EPI')
    testCondition = false;
end

UnitTest.tests(8).duration = toc(testTime);
UnitTest.tests(8).passed = testCondition;


%% Test run 9: Valid canonical input "EPI" + GE + 2D -> EPI unchanged, no side effects

UnitTest.tests(9).testname = 'Valid canonical EPI + GE + 2D -> EPI unchanged';

testTime = tic;

xQ = struct;
xQ.PulseSequenceType = 'EPI';
xQ.Vendor            = 'GE';
xQ.MRAcquisitionType = '2D';
xQ = xASL_adm_DefineASLReadout(xQ, false);

testCondition = true;
if ~strcmpi(xQ.PulseSequenceType, 'EPI')
    testCondition = false;
end
if ~strcmpi(xQ.Vendor, 'GE')
    testCondition = false;
end

UnitTest.tests(9).duration = toc(testTime);
UnitTest.tests(9).passed = testCondition;


%% Test run 10: Unknown PulseSequenceType (not fuzzy-recoverable) is stripped and inferred

UnitTest.tests(10).testname = 'Gibberish PulseSequenceType stripped, inferred from Vendor+MRAcq';

testTime = tic;

xQ = struct;
xQ.PulseSequenceType = 'zzzzzzzz'; % shares no letters in order with epi/grase/spiral
xQ.Vendor            = 'Siemens';
xQ.MRAcquisitionType = '3D';
xQ = xASL_adm_DefineASLReadout(xQ, false);

testCondition = true;
% Tier 2 fuzzy must reject this; field is stripped; Tier-3 inference sets GRASE
if ~isfield(xQ, 'PulseSequenceType') || ~strcmpi(xQ.PulseSequenceType, 'GRASE')
    testCondition = false;
end

UnitTest.tests(10).duration = toc(testTime);
UnitTest.tests(10).passed = testCondition;


%% Test run 11: Inference path — missing PulseSequenceType, 2D -> EPI

UnitTest.tests(11).testname = 'Missing PulseSequenceType + 2D -> inferred EPI';

testTime = tic;

xQ = struct;
xQ.Vendor            = 'Philips';
xQ.MRAcquisitionType = '2D';
xQ = xASL_adm_DefineASLReadout(xQ, false);

testCondition = true;
if ~isfield(xQ, 'PulseSequenceType') || ~strcmpi(xQ.PulseSequenceType, 'EPI')
    testCondition = false;
end

UnitTest.tests(11).duration = toc(testTime);
UnitTest.tests(11).passed = testCondition;


%% Test run 12: Inference path — missing PulseSequenceType, 3D + GE -> inferred spiral

UnitTest.tests(12).testname = 'Missing PulseSequenceType + 3D + GE -> inferred spiral';

testTime = tic;

xQ = struct;
xQ.Vendor            = 'GE';
xQ.MRAcquisitionType = '3D';
xQ = xASL_adm_DefineASLReadout(xQ, false);

testCondition = true;
if ~isfield(xQ, 'PulseSequenceType') || ~strcmpi(xQ.PulseSequenceType, 'spiral')
    testCondition = false;
end

UnitTest.tests(12).duration = toc(testTime);
UnitTest.tests(12).passed = testCondition;


%% Test run 13: Manufacturer substring match (e.g., "Siemens Healthineers" -> Siemens)

UnitTest.tests(13).testname = 'Manufacturer "Siemens Healthineers" substring-matches Vendor Siemens';

testTime = tic;

xQ = struct;
xQ.Manufacturer      = 'SIEMENS Healthineers';
xQ.MRAcquisitionType = '3D';
xQ = xASL_adm_DefineASLReadout(xQ, false);

testCondition = true;
if ~isfield(xQ, 'Vendor') || ~strcmpi(xQ.Vendor, 'Siemens')
    testCondition = false;
end

UnitTest.tests(13).duration = toc(testTime);
UnitTest.tests(13).passed = testCondition;


%% Test run 14: 2D + spiral triggers implausibility warning but does not erase the field

UnitTest.tests(14).testname = '2D + spiral implausible but field preserved';

testTime = tic;

xQ = struct;
xQ.PulseSequenceType = 'spiral'; % plausible value but implausible with 2D
xQ.Vendor            = 'GE';
xQ.MRAcquisitionType = '2D';
xQ = xASL_adm_DefineASLReadout(xQ, false);

testCondition = true;
% Consistency check warns only — must not erase
if ~isfield(xQ, 'PulseSequenceType') || ~strcmpi(xQ.PulseSequenceType, 'spiral')
    testCondition = false;
end

UnitTest.tests(14).duration = toc(testTime);
UnitTest.tests(14).passed = testCondition;


%% Restore warning state
warning(warnState);


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end