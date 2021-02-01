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
%               The idea is that this script can run independently from the
%               rest of ExploreASL, to enable unbiased and robust testing.
%
% EXAMPLE:      UnitTests = xASL_qc_UnitTesting;
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL

%% Initialize Paths

% Add testing directory
TestingPath = mfilename('fullpath');
filepath = fileparts(TestingPath);
addpath(filepath)

% Add all subfolders of xASL directory
xASL_dir = fileparts(filepath);
warning('off','all')
addpath(genpath(xASL_dir))
warning('on','all')

%% Clean Up
clc
clear

%% Get Test Repository
if usejava('desktop')
    TestRepository = uigetdir([],'Select test repository...');
else
    TestRepository = input('Insert test repository: ');
end

%% Test Workflow

% Unit test: xASL_tsvRead
UnitTests(1) = xASL_qc_UnitTest_function_tsvRead(TestRepository);

%% Print test results
clc
fprintf('================================= TEST RESULTS =================================\n')
for it = 1:numel(UnitTests)
    fprintf('TEST:\t\t%s\n',UnitTests(it).name)
    if UnitTests(it).passed
        fprintf('RESULTS:\t%s\n','Passed')
    else
        fprintf('RESULTS:\t%s\n','Failed')
    end
end
fprintf('================================================================================\n')







