function UnitTests = xASL_ut_UnitTesting
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
% 				This means you, as a developer, should not use/add ExploreASL
% 				functions within this script!
%
% EXAMPLE:      UnitTests = xASL_qc_UnitTesting;
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL

    %% Initialize Paths

    % Add testing directory
    TestingPath = mfilename('fullpath');    % Path of the current script
    filepath = fileparts(TestingPath);      % Testing directory
    xASL_dir = fileparts(filepath);         % All subfolders of xASL directory
    warning('off','all')                    % Turn of warnings (only add path related)
    addpath(genpath(xASL_dir))              % Add all scripts to path (including testing directory)
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
    if ~ischar(TestRepository)
        error('Test repository path required...');
    end

    %% Test Workflow

    % Unit test: xASL_tsvRead
    UnitTests(1) = xASL_ut_UnitTest_function_tsvRead(TestRepository);

    % Unit test: xASL_tsvWrite
    UnitTests(2) = xASL_ut_UnitTest_function_tsvWrite(TestRepository);

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



end



