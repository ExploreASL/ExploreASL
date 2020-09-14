function UnitTests = ExploreASL_Unit_Testing
%ExploreASL_Unit_Testing Main script to run all individual unit tests
%
% INPUT:        ...
%
% OUTPUT:       ...
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  ...
%
% EXAMPLE:      UnitTests = ExploreASL_Unit_Testing;
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2020 ExploreASL

%% Test Workflow
UnitTests = struct;

% Unit test A
[UnitTests(1).test, UnitTests(1).passed] = ExploreASL_Unit_Test_Template('Template Test', true);

% Unit test B
[UnitTests(2).test, UnitTests(2).passed] = ExploreASL_Unit_Test_Template('Template Test 2', false);




