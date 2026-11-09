function UnitTest = xASL_ut_function_xASL_bids_layout_IntendedFor(TestRepository)
%xASL_ut_function_xASL_bids_layout_IntendedFor Test M0 paths and BIDS URIs.
%
% INPUT:  TestRepository - Path to the test repository.
% OUTPUT: UnitTest       - Results for xASL_test_UnitTesting.
%
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE

warningState = warning('off', 'all');
warningCleanup = onCleanup(@() warning(warningState));
workingParent = fullfile(TestRepository, 'UnitTesting', 'working_directory');
xASL_adm_CreateDir(workingParent);
datasetRoot = tempname(workingParent);
xASL_adm_CreateDir(datasetRoot);
directoryCleanup = onCleanup(@() xASL_delete(datasetRoot, true));

perfDir = fullfile(datasetRoot, 'sub-001', 'ses-1', 'perf');
xASL_adm_CreateDir(perfDir);
xASL_io_WriteJson(fullfile(datasetRoot, 'dataset_description.json'), ...
    struct('Name', 'IntendedFor test', 'BIDSVersion', '1.8.0'));

for iRun = 1:2
    basename = sprintf('sub-001_ses-1_run-%02d', iRun);
    touch(fullfile(perfDir, [basename '_asl.nii.gz']));
    touch(fullfile(perfDir, [basename '_m0scan.nii.gz']));
    xASL_io_WriteJson(fullfile(perfDir, [basename '_asl.json']), struct('M0Type', 'Separate'));
    xASL_tsvWrite({'volume_type'; 'control'; 'label'}, fullfile(perfDir, [basename '_aslcontext.tsv']));
end

path1 = 'ses-1/perf/sub-001_ses-1_run-01_asl.nii.gz';
path2 = 'ses-1/perf/sub-001_ses-1_run-02_asl.nii.gz';
target1 = fullfile(datasetRoot, 'sub-001', path1);
target2 = fullfile(datasetRoot, 'sub-001', path2);
expectedFirst = {target1, target2};
firstM0 = fullfile(perfDir, 'sub-001_ses-1_run-01_m0scan.json');
secondM0 = fullfile(perfDir, 'sub-001_ses-1_run-02_m0scan.json');
xASL_io_WriteJson(firstM0, struct('IntendedFor', {{path1, path2}}));

% Name, second M0 sidecar, expected targets.
cases = {
    'String',               struct('IntendedFor', path2),                    {target2};
    'Single-element array', struct('IntendedFor', {{path2}}),                {target2};
    'Multiple targets',     struct('IntendedFor', {{path2, path1}}),         {target2, target1};
    'Mixed URI and path',   struct('IntendedFor', {{['bids::sub-001/' path2], path1}}), {target2, target1};
    'Session-relative',     struct('IntendedFor', 'perf/sub-001_ses-1_run-02_asl.nii.gz'), {target2};
    'Missing field',        struct(),                                      {};
    'Empty array',          struct('IntendedFor', {{}}),                     {};
    'Invalid type',         struct('IntendedFor', 42),                       {};
    'Mixed array',          struct('IntendedFor', {{path2, 42}}),            {};
    'Missing target',       struct('IntendedFor', 'ses-1/perf/missing.nii.gz'), {};
};

for iCase = 1:size(cases, 1)
    testTime = tic;
    xASL_io_WriteJson(secondM0, cases{iCase, 2});
    passed = false;
    try
        BIDS = bids.layout(datasetRoot);
        assert(numel(BIDS.subjects) == 1, 'Subject was not parsed.');
        perf = BIDS.subjects.perf;
        m0 = perf(strcmp({perf.type}, 'm0scan'));
        assert(numel(m0) == 2, 'Expected two M0 scans.');
        assert(isequal(m0(1).intended_for, expectedFirst), 'First M0 targets changed.');
        % The second M0 must not inherit targets from the first.
        assert(isequal(m0(2).intended_for, cases{iCase, 3}), 'Incorrect IntendedFor targets.');
        passed = true;
    catch ME
        fprintf('IntendedFor %s: %s\n', cases{iCase, 1}, ME.message);
    end
    UnitTest.tests(iCase).testname = cases{iCase, 1};
    UnitTest.tests(iCase).duration = toc(testTime);
    UnitTest.tests(iCase).passed = passed;
end

UnitTest = xASL_ut_CheckSubtests(UnitTest);
end

function touch(path)
% The parser checks filenames, not image contents.
    fid = fopen(path, 'w');
    if fid < 0
        error('Could not create test file: %s', path);
    end
    fclose(fid);
end
