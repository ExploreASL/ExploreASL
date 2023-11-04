function [summary, summaryTable] = xASL_test_Summarize(ResultsDir)
% xASL_test_Summarize Main script to combine all testing results into one summary
%
% FORMAT: [summary] = xASL_test_Summarize(ResultsDir)
%
% INPUT:        ResultsDir   - Pull up-to-date testing repository (BOOLEAN, REQUIRED)
%
% OUTPUT:       summary       - structure containing the summary of all test results        
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function is used after all test results have been generated in order to summarize a human readable file. 
%
% EXAMPLE:      [summary] = xASL_test_Summarize('/User/ASLTestResults/2022-20-01-results');
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2023 ExploreASL

MatFileList = xASL_adm_GetFileList(ResultsDir,'^*.mat$','List',[0 Inf], false);
TsvFileList = xASL_adm_GetFileList(ResultsDir,'^*.tsv$','List',[0 Inf], false);

summary  = struct;
tsvSummary = struct;

% Load in all mat files
for iList=1:length(MatFileList)
    % Get the full file path
    iFile = fullfile(ResultsDir, MatFileList{iList});
    % Filter out the name of the file without timestamp and extension.
    catName = extractBetween(MatFileList{iList}, 18, length(MatFileList{iList})-4);
    % Add the matfile to the summary if it's not another summary file.
    if (catName ~= "summary")
        summary(1).(catName{1}) = load(iFile, '*');
    end
end

% Load in all tsv files
for iList=1:length(TsvFileList)
    % Get files and directories
    iFile = fullfile(ResultsDir, TsvFileList{iList});
    iCell = xASL_tsvRead(iFile,false);
    % Filter out the name of the file without timestamp and extension.
    catName = extractBetween(TsvFileList{iList}, 18, length(TsvFileList{iList})-4);
    % Add the table to the summary 
    tsvSummary(1).(catName{1}) = iCell;
end

tsvFieldnames = fieldnames(tsvSummary);
summaryTable = struct;
for iList=1:length(tsvFieldnames)
    if strcmp(tsvFieldnames{iList}, 'unit_comparison')
        saveArray = tsvSummary.(tsvFieldnames{iList})(1,:);
        for iRow=2:size(tsvSummary.(tsvFieldnames{iList}),1)
            if ~tsvSummary.(tsvFieldnames{iList}){iRow,3}
                saveRow = tsvSummary.(tsvFieldnames{iList})(iRow,:);
                saveArray = [saveArray; saveRow];
            end
        end
        summaryTable(1).(tsvFieldnames{iList}) = saveArray;
    elseif strcmp(tsvFieldnames{iList},'flavor_comparison')
        saveArray = tsvSummary.(tsvFieldnames{iList})(1,1:4);
        for iRow=2:size(tsvSummary.(tsvFieldnames{iList}),1)
            if true % ~tsvSummary.(tsvFieldnames{iList}){iRow,5}% FUNCTIONALITY NOT ACTIVE
                saveRow = tsvSummary.(tsvFieldnames{iList})(iRow,1:4);
                saveArray = [saveArray; saveRow];
            end
        end
        summaryTable(1).(tsvFieldnames{iList}) = saveArray;
    elseif regexp(tsvFieldnames{iList}, '^ASL_.*ComparisonTable.*$')
        saveArray = tsvSummary.(tsvFieldnames{iList})(1,:);
        for iRow=2:size(tsvSummary.(tsvFieldnames{iList}),1)
            if sum([tsvSummary.(tsvFieldnames{iList}){iRow,2:end}]) < 10
                saveRow = tsvSummary.(tsvFieldnames{iList})(iRow,:);
                saveArray = [saveArray; saveRow];
            end
        end
        summaryTable(1).(tsvFieldnames{iList}) = saveArray;
    end
end

% Save path
TimeString = datestr(now,'yyyy-mm-dd_HH_MM');
TimeString(end-2) = 'h';
savePathMat = fullfile(ResultsDir, [TimeString, '_summary.mat']);
savePathTSV = fullfile(ResultsDir, [TimeString, '_summary.tsv']);

% Save files
%xASL_tsvWrite(UnitTestsCells,savePathTSV);
save(savePathMat,'summary');



%% Print test results
%clc
fprintf('====================================== Summary Mat ======================================\n\n')
disp(summary);
fprintf('====================================== Summary TSV ======================================\n\n')
disp(tsvSummary);
fprintf('====================================== Fails Only  ======================================\n\n')
disp(summaryTable);
fprintf('==========================================================================================\n')

end

