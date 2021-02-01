function UnitTests = xASL_qc_UnitTesting
%xASL_qc_UnitTesting Main script to run all individual unit tests
%
% INPUT:        n/a
%
% OUTPUT:       UnitTests structure
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function is used to run all individual unit tests. To
%               define a unit test, please use the xASL_qc_UnitTest_Template.
%
% EXAMPLE:      UnitTests = xASL_qc_UnitTesting;
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2020 ExploreASL

%% Test Workflow

% Unit test A
UnitTests(1) = xASL_qc_UnitTest_Template;

% Unit test B
UnitTests(2) = xASL_qc_UnitTest_Template;




