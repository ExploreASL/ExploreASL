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
[UnitTests(1).name,UnitTests(1).module,UnitTests(1).submodule,UnitTests(1).passed,UnitTests(1).tests] = ExploreASL_Unit_Test_Template;

% Unit test B
[UnitTests(2).name,UnitTests(2).module,UnitTests(2).submodule,UnitTests(2).passed,UnitTests(2).tests] = ExploreASL_Unit_Test_Template;




