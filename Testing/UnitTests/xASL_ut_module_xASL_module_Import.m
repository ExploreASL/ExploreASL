function UnitTest = xASL_ut_module_xASL_module_Import(TestRepository)
%xASL_ut_module_xASL_module_Import Individual unit test for xASL_module_Import
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
% EXAMPLE:      UnitTests(1) = xASL_ut_module_xASL_module_Import(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________



%% Test run 1

UnitTest.tests(1).testname = 'Failed NII2BIDS does not mark import ready';
testTime = tic;

% Prepare a minimal NII2BIDS input that fails while reading the ASL JSON
pathTestPatient = fullfile(TestRepository,'UnitTesting','working_directory','ImportFailedStatus');
pathTempASL = fullfile(pathTestPatient,'derivatives','ExploreASL','temp','Sub1','ASL_1');
pathImportLock = fullfile(pathTestPatient,'derivatives','ExploreASL','lock','xASL_module_Import','Sub1','xASL_module_Import');
xASL_adm_CreateDir(fullfile(pathTestPatient,'sourcedata'));
xASL_adm_CreateDir(pathTempASL);
xASL_adm_CreateDir(pathImportLock);

fid = fopen(fullfile(pathTempASL,'ASL4D.json'),'wt');
fprintf(fid,'{ this is not valid json');
fclose(fid);
% Simulate a successful DCM2NII
fid = fopen(fullfile(pathImportLock,'010_DCM2NII.status'),'wt');
fclose(fid);
% Deliberate stale ready status
fid = fopen(fullfile(pathImportLock,'999_ready.status'),'wt');
fclose(fid);

ExploreASL(pathTestPatient,[0 1 0],0,0,1,1);

% Input should remain untouched and completion states must not be written
testCondition = exist(fullfile(pathTempASL,'ASL4D.json'),'file') == 2 && ...
    exist(fullfile(pathImportLock,'010_DCM2NII.status'),'file') == 2 && ...
    exist(fullfile(pathImportLock,'020_NII2BIDS.status'),'file') ~= 2 && ...
    exist(fullfile(pathImportLock,'999_ready.status'),'file') ~= 2;

% Clean up
xASL_delete(pathTestPatient,true)
UnitTest.tests(1).duration = toc(testTime);
UnitTest.tests(1).passed = testCondition;


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end
