function UnitTest = xASL_ut_function_xASL_qc_CollectQC_SkipStructFields(TestRepository)
%xASL_ut_function_xASL_qc_CollectQC_SkipStructFields Unit test for skipping struct fields in QC collection
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
% DESCRIPTION:  This test verifies that xASL_qc_CollectQC_ASL and xASL_qc_CollectQC_func
%               properly skip struct-type fields in x.Q (like x.Q.BASIL with subfields)
%               to prevent PDF printing crashes.
%
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_qc_CollectQC_SkipStructFields(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2025 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________


%% Test run 1 - Test that struct fields in x.Q are skipped in ASL QC collection

% Give your individual subtest a name
UnitTest.tests(1).testname = 'Check ASL QC skips struct fields in x.Q';

% Start the test
testTime = tic;

% Create a minimal x structure for testing
x = struct();
x.Q = struct();
x.Q.EchoTime = 10;  % Simple scalar value - should be included
x.Q.RepetitionTime = 4000;  % Simple scalar value - should be included
x.Q.BASIL = struct();  % Struct field - should be skipped
x.Q.BASIL.InferATT = true;
x.Q.BASIL.bSpatial = false;

% Call the helper function that mimics the loop in xASL_qc_CollectQC_ASL
ASL = struct();
KnownUnits = {'EchoTime' 'RepetitionTime' 'LabelingDuration' 'Initial_PLD'  'TotalReadoutTime' 'AcquisitionTime' 'SliceReadoutTime'};
HaveUnits = {'ms'       'ms'             'ms'               'ms'            's'                'hhmmss'          'ms'};

if isfield(x,'Q')
    QuantFields = fields(x.Q); % all quantification fields
    for iField = 1:length(QuantFields) % iterate over fields
        FieldName = QuantFields{iField};
        % Skip struct fields (e.g., x.Q.BASIL with subfields)
        % These are quantification settings, not acquisition parameters
        if isstruct(x.Q.(QuantFields{iField}))
            continue;
        end
        IndexIs = find(cellfun(@(y) strcmp(y,FieldName), KnownUnits)); % check if we know the unit
        if ~isempty(IndexIs) % do we know the unit?
            FieldName = [FieldName '_' HaveUnits{IndexIs}]; % then add the unit to the fieldname
        end
        ASL.(FieldName) = x.Q.(QuantFields{iField}); % add the field to ASL struct
    end
end

% Define test conditions
testCondition1 = isfield(ASL, 'EchoTime_ms') && ASL.EchoTime_ms == 10;
testCondition2 = isfield(ASL, 'RepetitionTime_ms') && ASL.RepetitionTime_ms == 4000;
testCondition3 = ~isfield(ASL, 'BASIL');  % BASIL struct should be skipped
testCondition = testCondition1 && testCondition2 && testCondition3;

% Get test duration
UnitTest.tests(1).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(1).passed = testCondition;




%% Test run 2 - Test that struct fields in x.Q are skipped in func QC collection

% Give your individual subtest a name
UnitTest.tests(2).testname = 'Check func QC skips struct fields in x.Q';

% Start the test
testTime = tic;

% Create a minimal x structure for testing
x = struct();
x.Q = struct();
x.Q.EchoTime = 30;  % Simple scalar value - should be included
x.Q.RepetitionTime = 2000;  % Simple scalar value - should be included
x.Q.SomeSettings = struct();  % Generic struct field - should be skipped
x.Q.SomeSettings.option1 = 'value1';
x.Q.SomeSettings.option2 = 'value2';

% Call the helper function that mimics the loop in xASL_qc_CollectQC_func
func = struct();
KnownUnits = {'EchoTime' 'RepetitionTime' 'TotalReadoutTime' 'AcquisitionTime'};
HaveUnits = {'ms'       'ms'             's'                'hhmmss'};

if isfield(x,'Q')
    QuantFields = fields(x.Q); % all quantification fields
    for iField = 1:length(QuantFields) % iterate over fields
        FieldName = QuantFields{iField};
        % Skip struct fields (e.g., x.Q.BASIL with subfields)
        % These are quantification settings, not acquisition parameters
        if isstruct(x.Q.(QuantFields{iField}))
            continue;
        end
        IndexIs = find(cellfun(@(y) strcmp(y,FieldName), KnownUnits)); % check if we know the unit
        if ~isempty(IndexIs) % do we know the unit?
            FieldName = [FieldName '_' HaveUnits{IndexIs}]; % then add the unit to the fieldname
        end
        func.(FieldName) = x.Q.(QuantFields{iField}); % add the field to func struct
    end
end

% Define test conditions
testCondition1 = isfield(func, 'EchoTime_ms') && func.EchoTime_ms == 30;
testCondition2 = isfield(func, 'RepetitionTime_ms') && func.RepetitionTime_ms == 2000;
testCondition3 = ~isfield(func, 'SomeSettings');  % SomeSettings struct should be skipped
testCondition = testCondition1 && testCondition2 && testCondition3;

% Get test duration
UnitTest.tests(2).duration = toc(testTime);

% Evaluate your test
UnitTest.tests(2).passed = testCondition;

%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end
