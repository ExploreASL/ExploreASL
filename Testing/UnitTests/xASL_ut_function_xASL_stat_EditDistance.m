function UnitTest = xASL_ut_function_xASL_stat_EditDistance(TestRepository)
%xASL_ut_function_xASL_stat_EditDistance Unit test for xASL_stat_EditDistance
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
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_stat_EditDistance(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________



%% ---- Group A: Basic distance properties --------------------------------

% Subtest 1: Zero distance for equal strings
UnitTest.tests(1).testname = 'Zero distance for identical strings';
testTime = tic;
d = xASL_stat_EditDistance('grase', 'grase');
UnitTest.tests(1).passed = (d == 0);
UnitTest.tests(1).duration = toc(testTime);

% Subtest 2: Distance equals the length of the non-empty string when the other is empty
UnitTest.tests(2).testname = 'Empty vs non-empty: distance equals string length';
testTime = tic;
d = xASL_stat_EditDistance('', 'epi');
UnitTest.tests(2).passed = (d == 3);
UnitTest.tests(2).duration = toc(testTime);

% Subtest 3: Both strings empty yields zero
UnitTest.tests(3).testname = 'Both empty: distance is zero';
testTime = tic;
d = xASL_stat_EditDistance('', '');
UnitTest.tests(3).passed = (d == 0);
UnitTest.tests(3).duration = toc(testTime);


%% ---- Group B: Single-operation correctness -----------------------------

% Subtest 4: One character changed (substitution only)
UnitTest.tests(4).testname = 'Single character substitution';
testTime = tic;
d = xASL_stat_EditDistance('epi', 'epo');
UnitTest.tests(4).passed = (d == 1);
UnitTest.tests(4).duration = toc(testTime);

% Subtest 5: One character appended (insertion only)
UnitTest.tests(5).testname = 'Single character insertion at end';
testTime = tic;
d = xASL_stat_EditDistance('grase', 'grasee');
UnitTest.tests(5).passed = (d == 1);
UnitTest.tests(5).duration = toc(testTime);

% Subtest 6: One character removed (deletion only)
UnitTest.tests(6).testname = 'Single character deletion';
testTime = tic;
d = xASL_stat_EditDistance('spiral', 'spial');
UnitTest.tests(6).passed = (d == 1);
UnitTest.tests(6).duration = toc(testTime);

% Subtest 7: Two adjacent characters swapped (transposition only)
UnitTest.tests(7).testname = 'Adjacent-character transposition: ri -> ir';
testTime = tic;
d = xASL_stat_EditDistance('spiral', 'sprial');
UnitTest.tests(7).passed = (d == 1);
UnitTest.tests(7).duration = toc(testTime);


%% ---- Group C: Damerau vs Levenshtein distinction -----------------------

% Subtest 8: Transposition distance is 1 (Damerau), not 2 (Levenshtein)
% 'Seimens' and 'Siemens' differ only by adjacent swap of 'ei' -> 'ie'
UnitTest.tests(8).testname = 'Manufacturer typo: Seimens to Siemens is transposition, distance 1';
testTime = tic;
d = xASL_stat_EditDistance('Seimens', 'Siemens');
UnitTest.tests(8).passed = (d == 1);
UnitTest.tests(8).duration = toc(testTime);

% Subtest 9: Non-adjacent swap is NOT a single transposition
% 'abcdef' vs 'abcfde': swap d<->f would need d->e (subst) + e->f (subst) + remove extra = 2 ops
% Actually: abcdef -> abcfde: at pos4-5-6, 'def' vs 'fde'. 
%   delete 'd' (1), now 'abcef' vs 'abcfde', 
%   'abcef' -> 'abcfe' (transpose e,f: 1), 'abcfe' -> 'abcfde' (insert d: 1) = 3
%   Or: delete 'e' at pos5 (1), 'abcdf' vs 'abcfde', insert 'e' at pos5 (1), = 2
%   Or: substitute d->f (1), e->d (1) = 2
%   Best = 2.
UnitTest.tests(9).testname = 'Non-adjacent characters differ: distance > 1';
testTime = tic;
d = xASL_stat_EditDistance('abcdef', 'abcfde');
UnitTest.tests(9).passed = (d == 2);
UnitTest.tests(9).duration = toc(testTime);


%% ---- Group D: Longer / multi-operation sequences -----------------------

% Subtest 10: Standard example from literature
UnitTest.tests(10).testname = 'Literature example: kitten to sitting is distance 3';
testTime = tic;
d = xASL_stat_EditDistance('kitten', 'sitting');
UnitTest.tests(10).passed = (d == 3);
UnitTest.tests(10).duration = toc(testTime);

% Subtest 11: Completely different strings (no shared characters)
UnitTest.tests(11).testname = 'Completely different strings: abc vs xyz is distance 3';
testTime = tic;
d = xASL_stat_EditDistance('abc', 'xyz');
UnitTest.tests(11).passed = (d == 3);
UnitTest.tests(11).duration = toc(testTime);

% Subtest 12: Prefix relationship (one string is a prefix of the other)
UnitTest.tests(12).testname = 'Prefix relationship: ep2d to epi requires 2 edits';
testTime = tic;
d = xASL_stat_EditDistance('ep2d', 'epi');
UnitTest.tests(12).passed = (d == 2);
UnitTest.tests(12).duration = toc(testTime);


%% ---- Group E: Motivating real-world ASL use case -----------------------

% Subtest 13: PulseSequenceType typo with dimension prefix stripped
% '3D_SPRIAL' lowercased to '3d_sprial', stripped of '3d_' gives 'sprial',
% then fuzzy match to 'spiral' with distance 1
UnitTest.tests(13).testname = 'ASL real case: 3D_SPRIAL normalized matches spiral (distance 1 from sprial)';
testTime = tic;
pstLower   = lower('3D_SPRIAL');
pstStripped = regexprep(pstLower, '^(\d+d[_\-]?)', '');
d = xASL_stat_EditDistance(pstStripped, 'spiral');
UnitTest.tests(13).passed = (d == 1);
UnitTest.tests(13).duration = toc(testTime);

% Subtest 14: Another common typo — spiril is distance 1 from spiral
UnitTest.tests(14).testname = 'ASL real case: spiril (letter swap) is distance 1 from spiral';
testTime = tic;
d = xASL_stat_EditDistance('spiril', 'spiral');
UnitTest.tests(14).passed = (d == 1);
UnitTest.tests(14).duration = toc(testTime);


%% ---- Group F: Symmetry and cost parameter tests ------------------------

% Subtest 15: Symmetry — d(s1,s2) == d(s2,s1) with uniform costs
UnitTest.tests(15).testname = 'Symmetry property with uniform costs';
testTime = tic;
d12 = xASL_stat_EditDistance('spiral', 'epiral');
d21 = xASL_stat_EditDistance('epiral', 'spiral');
UnitTest.tests(15).passed = (d12 == d21) && (d12 == 1);
UnitTest.tests(15).duration = toc(testTime);

% Subtest 16: Custom costs — high replacement cost forces del+ins instead
UnitTest.tests(16).testname = 'High replacement cost favors deletion + insertion';
testTime = tic;
% 'a' -> 'b': with ReplCost=10, cheaper to delete 'a' (1) + insert 'b' (1) = 2
d = xASL_stat_EditDistance('a', 'b', 1, 1, 10);
UnitTest.tests(16).passed = (d == 2);
UnitTest.tests(16).duration = toc(testTime);


%% End of testing
UnitTest = xASL_ut_CheckSubtests(UnitTest);

end