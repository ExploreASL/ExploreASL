function RESULTS = xASL_qc_UnitTesting
%xASL_qc_UnitTesting Script to run all module and submodule tests
%
% FORMAT:       RESULTS = xASL_qc_UnitTesting
% 
% INPUT:        None
%
% OUTPUT:       RESULTS structure
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Script to run all module and submodule tests. Please start
%               this script from the ExploreASL directory and select a test
%               directory. A temporary TestFolder is created. If it already
%               exists it will be removed.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLES:     RESULTS = xASL_qc_UnitTesting;
% __________________________________
% Copyright 2015-2020 ExploreASL

% Improve command window output
BreakString = [repmat('=',1,100),'\n'];

%% Update GIT
try
    fprintf('GIT: ');
    system('git fetch','-echo');
    system('git pull','-echo');
catch
    warning('It seems that your working directory is not the ExploreASL directory...');
end


%% GET DIRECTORIES
test_parameter_file = fullfile(pwd,'Development','ExploreASL_UnitTesting','xASL_test_parameters.json');
if xASL_exist(test_parameter_file)
    val = jsondecode(fileread(test_parameter_file));
else
    % If the file does not exist yet, it needs to be created
    val = struct;
end

% Check if the current directory is the ASL directory
if xASL_exist(fullfile(pwd,'ExploreASL_Master'))
    val.xASLdir = 'pwd';
else
    val.xASLdir = uigetdir(pwd, 'Select ExploreASL directory...');
end

% Define test directory
val.testDir = uigetdir(pwd, 'Select testing directory...');

% Write the directories to the JSON file
valNew = jsonencode(val);
fid = fopen(test_parameter_file, 'w');
fprintf(fid, '%s', valNew);
fclose(fid);

%% START TESTING
fprintf(BreakString);


%% RUN TESTS: INITIALIZATION
RESULTS.INIT = runtests('xASL_Initialize_test'); fprintf(BreakString);


%% RUN TESTS: MODULES
% RESULTS.MODULES.STRUCTURAL = runtests('xASL_module_Structural_test'); fprintf(BreakString);
% RESULTS.MODULES.ASL = runtests('xASL_module_ASL_test'); fprintf(BreakString);
% RESULTS.MODULES.POPULATION = runtests('xASL_module_Population_test'); fprintf(BreakString);


%% RUN TESTS: SUBMODULES
RESULTS.SUBMODULES.xASL_wrp_LinearReg_T1w2MNI = runtests('xASL_wrp_LinearReg_T1w2MNI_test'); fprintf(BreakString);



%% RESULTS PDF

% Get path
resultPath = fullfile(pwd,'Development','ExploreASL_UnitTesting','Results_QC_UnitTesting.pdf');
fg = spm_figure('Create','Graphics','visible','off');
set(fg,'windowstyle','normal');
set(fg,'units','normalized');
axDoc = axes('Position',[0 0 1 1]);
set(axDoc, 'visible', 'off')
fontsize = 8;

% Display results
text(0,0.99, 'RESULTS: QC UNIT TESTING', 'FontSize', fontsize+4, 'FontWeight', 'bold', 'Interpreter', 'none', 'Parent', axDoc);
text(0,0.96, 'INITIALIZATION', 'FontSize', fontsize+2, 'FontWeight', 'bold', 'Interpreter', 'none', 'Parent', axDoc);
printTestResults(RESULTS.INIT,0.00,0.94,axDoc)      % Print individual test results

text(0,0.70, 'SUBMODULES', 'FontSize', fontsize+2, 'FontWeight', 'bold', 'Interpreter', 'none', 'Parent', axDoc);
printTestResults(RESULTS.SUBMODULES.xASL_wrp_LinearReg_T1w2MNI,0.00,0.68,axDoc)      % Print individual test results

% Print PDF
print(fg, '-dpdf', '-r600', resultPath);


end

% Function to print the results of a single test
function printTestResults(testStruct,startX,startY,curAxes)

    % Configuration
    diffX = 0.20;
    diffY = 0.02;
    fontsize = 8;
    tmpDiff = 0;
    
    % Iterate over multiple tests
    for i=1:length(testStruct)
        text(startX, startY-tmpDiff, 'Name', 'FontSize', fontsize,...
             'Interpreter', 'none', 'Parent', curAxes);
        text(startX+diffX, startY-tmpDiff, string(testStruct(1,i).Name), 'FontSize', fontsize,...
             'Interpreter', 'none', 'Parent', curAxes);

        text(startX, startY-diffY-tmpDiff, 'Passed', 'FontSize', fontsize,...
             'Interpreter', 'none', 'Parent', curAxes);
        text(startX+diffX, startY-diffY-tmpDiff, string(testStruct(1,i).Passed), 'FontSize', fontsize,...
             'Interpreter', 'none', 'Parent', curAxes);

        text(startX, startY-diffY*2-tmpDiff, 'Duration', 'FontSize', fontsize,...
             'Interpreter', 'none', 'Parent', curAxes);
        text(startX+diffX, startY-diffY*2-tmpDiff, strcat(string(testStruct(1,i).Duration),'s'), 'FontSize', fontsize,...
             'Interpreter', 'none', 'Parent', curAxes);
         
        tmpDiff = tmpDiff+0.06;
    end


end
