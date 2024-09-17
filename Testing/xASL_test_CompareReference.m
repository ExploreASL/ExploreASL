function [result] = xASL_test_CompareReference(pathReference, pathResults, pathDestination)
%xASL_test_CompareReference compare the output of ExploreASL TestDataSets to reference values
%
% FORMAT: [result] = xASL_test_CompareReference(pathReference, pathResults [, pathDestination])
% 
% INPUT:
%   pathReference - path to file containing all reference values to compare to (REQUIRED)
%   pathResults   - path to folder where the test results are stored (REQUIRED)
%   pathDestination  - path to folder where the test results moved to (OPTIONAL)
%                      Defaults to saving in pathResults. 
%                 
% OUTPUT:
%   result - Table containing all results from the test runs
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function will take the output from a processed version of TestDataSets and make a set of tables comparing it to reference values.
%              The reference values are strored in ExploreASL/Testing/ReferenceValues.tsv
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:
%
%      [result] = xASL_test_CompareReference('<dir>/ExploreASL/Testing//ReferenceValues.tsv', '<dir>/TestDataSetsTemp', '<dir>/TestResults');
% __________________________________
% Copyright (c) 2015-2023 ExploreASL 
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

% Validate the input options
if nargin < 2 || isempty(pathReference) || isempty(pathResults)
    error('Require three input parameters pathReference, pathResults');
end
if nargin < 3 || isempty(pathDestination)
    pathDestination = pathResults;
end
result = struct;
Dlist = xASL_adm_GetFileList(pathResults,'^.*$','List',[0 Inf], true);
[result.ResultsTable, result.ResultTableFile,  result.SaveFile] = xASL_test_DetermineResultsTable(pathResults, Dlist, pathDestination);
[result.ReferenceTables , result.ReferenceTable] = xASL_qc_LoadRefTable(pathReference);
[result.ResultsComparison, result.ResultsDifference] = xASL_qc_CompareTables(pathResults, result.ReferenceTable, result.ResultsTable, pathDestination);  
end    
    
% Load Reference Table
function [ReferenceTables,ReferenceTable] = xASL_qc_LoadRefTable(pathRefTable)
    % Load TSV file
    ReferenceTables = xASL_tsvRead(pathRefTable);
    
    iRow = 1;
    while iRow<=size(ReferenceTables,1)
        if ~isempty(regexp(ReferenceTables{iRow,1},'^xASL_', 'once'))
            versionXASL = ReferenceTables{iRow,1};
            operatingSystem = ReferenceTables{iRow+1,1};
            fprintf('Version: %s\n', versionXASL);
            fprintf('OS:      %s\n', operatingSystem);
            ReferenceTables(iRow,:) = []; % Remove row 1
            ReferenceTables(iRow,:) = []; % Remove row 2
            ReferenceTable.(versionXASL) = ReferenceTables(iRow:iRow+10,:);
        end
        iRow=iRow+1;
    end
end
%% Compare Results and Reference Tables
function [ResultsComparison, ResultsDifference] = xASL_qc_CompareTables(TestDir, ReferenceTable, ResultsTable, SaveDir)
    % Compare tables (skip first row)
    ResultsComparison = ReferenceTable;
    ResultsDifference = ReferenceTable;
    
    % Define Results comparison name
    TimeString = datestr(now,'yyyy-mm-dd_HH_MM');
    TimeString(end-2) = 'h';
    
    % Iterate over versions
    versionsXASL = fieldnames(ResultsComparison);
    for iVersion = 1:size(versionsXASL,1)
        % Get current reference table
        currentReferenceTable = ReferenceTable.(versionsXASL{iVersion});
        for iRow = 2:size(currentReferenceTable,1)
            ResultsComparison.(versionsXASL{iVersion})(iRow,2:end) = {NaN};
            ResultsDifference.(versionsXASL{iVersion})(iRow,2:end) = {NaN};
            if strcmp(currentReferenceTable{iRow,1},ResultsTable{iRow,1})
                % Iterate over results
                for iColumn = 2:size(currentReferenceTable,2)
                    refValue = xASL_str2num(currentReferenceTable{iRow,iColumn});
                    resValue = xASL_str2num(ResultsTable{iRow,iColumn});
					if isnan(refValue) && isnan(resValue)
						% NaNs are assigned to strings, if both are strings or NaNs then the comparison is fine
						ResultsComparison.(versionsXASL{iVersion})(iRow,iColumn) = {1};
                        ResultsDifference.(versionsXASL{iVersion})(iRow,iColumn) = {'NaN'};
					else
						% Check if difference is smaller than 0.1% of reference value
						if abs(refValue-resValue)<(abs(refValue)*0.01)
							ResultsComparison.(versionsXASL{iVersion})(iRow,iColumn) = {1};
                            ResultsDifference.(versionsXASL{iVersion})(iRow,iColumn) = {[xASL_num2str(100*abs(refValue-resValue)/(abs(refValue)), '%.3f') '%']};
						else
							ResultsComparison.(versionsXASL{iVersion})(iRow,iColumn) = {0};
                            ResultsDifference.(versionsXASL{iVersion})(iRow,iColumn) = {[xASL_num2str(100*abs(refValue-resValue)/(abs(refValue)), '%.3f') '%']};
						end
					end
                end
            else
                fprintf('Dataset names do not match...\n');
                fprintf('Reference: %s\n', currentReferenceTable{iRow,1});
                fprintf('Results:   %s\n', ResultsTable{iRow,1});
            end
        end
    end
    
    % Check: single value passed/failed
    versionsXASL = fieldnames(ResultsComparison);
    for iVersion = 1:size(versionsXASL,1)
        % Get current reference table
        currentResults = ResultsComparison.(versionsXASL{iVersion});
        resultValues = currentResults(2:end,2:end);
        % Check if there are "failed" values
        passed = true;
        for iRow=1:size(resultValues,1)
            for iColumn=1:size(resultValues,2)
                %fprintf('%d\n', resultValues{iRow,iColumn});
                if resultValues{iRow,iColumn}==0
                    passed = false; % At least one value is "failed"
                end
            end
        end
        % Store value
        ResultsComparison.([versionsXASL{iVersion} '_Passed']) = passed;
        if passed
            fprintf('Version: %s\nPassed:  %s\n', versionsXASL{iVersion}, 'true');
            ComparisonTableTsv = [TimeString, versionsXASL{iVersion}, '_ComparisonTable_Passed.tsv'];
        else
            fprintf('Version: %s\nPassed:  %s\n', versionsXASL{iVersion}, 'false');
            ComparisonTableTsv = [TimeString, versionsXASL{iVersion}, '_ComparisonTable_Failed.tsv'];
        end
        
        DifferencesTableTsv = [TimeString, versionsXASL{iVersion}, '_DifferenceTable.tsv'];
        xASL_tsvWrite(ResultsComparison.(versionsXASL{iVersion}), fullfile(SaveDir, ComparisonTableTsv));
        xASL_tsvWrite(ResultsDifference.(versionsXASL{iVersion}), fullfile(SaveDir, DifferencesTableTsv));
    end
end
%% Determine the results table
function [ResultsTable,ResultTableFile,SaveFile] = xASL_test_DetermineResultsTable(TestDir, Dlist, SaveDir)
    % Define results table name & fields
    ResultTableName = datestr(now,'yyyy-mm-dd_HH_MM');
    ResultTableName(end-2) = 'h';
    ResultsTable = {'Data', 'mean_qCBF_TotalGM' 'median_qCBF_TotalGM' 'median_qCBF_DeepWM' 'CoV_qCBF_TotalGM' 'GMvol' 'WMvol' 'CSFvol' 'PipelineCompleted' 'TC_ASL_Registration' 'TC_M0_Registration'};
    
    % Initialize again to make sure that the TemplateDir exists
    x = ExploreASL;
    
    fprintf('Reading & parsing results:   ');
    for iList=1:length(Dlist) % iterate over example datasets
        xASL_TrackProgress(iList, length(Dlist));
        ResultsTable{1+iList,1} = Dlist{iList};
        bidsDir = fullfile(TestDir, Dlist{iList});
        ExploreASLDir = fullfile(bidsDir, 'derivatives', 'ExploreASL');
        PopulationDir = fullfile(ExploreASLDir, 'Population');
        StatsDir = fullfile(PopulationDir, 'Stats');
        VolumeDir = fullfile(PopulationDir, 'TissueVolume');
        
        clear ResultsFile
        ResultFile{1} = xASL_adm_GetFileList(StatsDir,'(?i)^mean_qCBF.*TotalGM.*PVC2\.tsv$','FPList');
        ResultFile{2} = xASL_adm_GetFileList(StatsDir,'(?i)^median_qCBF.*TotalGM.*PVC0\.tsv$','FPList');
        ResultFile{3} = xASL_adm_GetFileList(StatsDir,'(?i)^median_qCBF.*DeepWM.*PVC0\.tsv$','FPList');
        ResultFile{4} = xASL_adm_GetFileList(StatsDir,'(?i)^CoV_qCBF.*TotalGM.*PVC0\.tsv$','FPList');
        ResultFile{5} = xASL_adm_GetFileList(VolumeDir,'(?i)^TissueVolume.*\.tsv$','FPList');
        
        for iFile=1:length(ResultFile) % iterate over ROI results
            % Make sure one individual file does not crash the table generation
            try
                if length(ResultFile{iFile})<1
                    ResultsTable{1+iList,1+iFile} = 'empty';
                    ResultsTable{1+iList,2+iFile} = 'empty';
                    ResultsTable{1+iList,3+iFile} = 'empty';
                elseif iFile<5 % check the ASL parameters
                    [~, TempTable] = xASL_bids_csv2tsvReadWrite(ResultFile{iFile}{end});
                    ResultsTable{1+iList,1+iFile} = TempTable{3,end-2};
                else % check the volumetric parameters
                    [~, TempTable] = xASL_bids_csv2tsvReadWrite(ResultFile{iFile}{end});
                    % Backward compatibility:
                    % Volumetrics used to be saved as '_(L)', but this is converted by xASL_io_ReadJson to its HEX counterpart '_0x28L0x29'
                    % Now we always save to _L to avoid this. For backward compatibility we still check the old options here
                    IndexGM = find(cellfun(@(y) ~isempty(regexpi(y,'(GM_volume_L|GM_volume_(L)|GM_volume_0x28L0x29)')), TempTable(1,:)));
                    IndexWM = find(cellfun(@(y) ~isempty(regexpi(y,'(WM_volume_L|WM_volume_(L)|WM_volume_0x28L0x29)')), TempTable(1,:)));
                    IndexCSF = find(cellfun(@(y) ~isempty(regexpi(y,'(CSF_volume_L|CSF_volume_(L)|CSF_volume_0x28L0x29)')), TempTable(1,:)));
                    if ~isempty(IndexGM)
                        ResultsTable{1+iList,1+iFile} = TempTable{2, IndexGM};
                    else
                        ResultsTable{1+iList,1+iFile} = 'n/a';
                    end
                    if ~isempty(IndexWM)
                        ResultsTable{1+iList,2+iFile} = TempTable{2, IndexWM};
                    else
                        ResultsTable{1+iList,2+iFile} = 'n/a';
                    end
                    if ~isempty(IndexCSF)
                        ResultsTable{1+iList,3+iFile} = TempTable{2, IndexCSF};
                    else
                        ResultsTable{1+iList,3+iFile} = 'n/a';
                    end
                end
            catch ME
                % Something went wrong, we set all values to n/a
                fprintf('%s\n', ME.message);
                ResultsTable{1+iList,1+iFile} = 'n/a';
                ResultsTable{1+iList,2+iFile} = 'n/a';
                ResultsTable{1+iList,3+iFile} = 'n/a';
            end
        end
        % check if there are missing lock files
        if exist(fullfile(bidsDir,'Missing_Lock_files.csv'),'file')
            % pipeline not completed
            ResultsTable{1+iList,4+length(ResultFile)} = 0;
        else
            % pipeline completed
            ResultsTable{1+iList,4+length(ResultFile)} = 1;
        end
        
        % Get registration performance
        ResultsTable{1+iList,5+length(ResultFile)} = 'n/a';
        ResultsTable{1+iList,6+length(ResultFile)} = 'n/a';
        
        qcFile = xASL_adm_GetFileList(ExploreASLDir, '^QC_collection.*\.json$', 'FPListRec');
        if ~isempty(qcFile)
            json = xASL_io_ReadJson(qcFile{1});
            if isfield(json, 'ASL')
                RunsAre = fields(json.ASL);
                bCBFFound = false;
                bM0Found = false;
                for iRun=1:length(RunsAre)
                    % now we try to find these TC parameters for any of the runs,
                    % we want this to be the first for backward compatibility
                    if ~bCBFFound && isfield(json.ASL.(RunsAre{iRun}), 'TC_CBF2template')
                        ResultsTable{1+iList,5+length(ResultFile)} = json.ASL.(RunsAre{iRun}).TC_CBF2template;
                        bCBFFound = true;
                    end
                    if ~bM0Found && isfield(json.ASL.(RunsAre{iRun}), 'TC_M02template')
                        ResultsTable{1+iList,6+length(ResultFile)} = json.ASL.(RunsAre{iRun}).TC_M02template;
                        bM0Found = true;
                    end
                end
            end
        end
    end
    fprintf('\n');
    
    % Save results
    ResultTableFile = [ResultTableName,'_ResultsTable.mat'];
    ResultTableTsv = [ResultTableName,'_ResultsTable.tsv'];
    SaveFile = fullfile(SaveDir, ResultTableFile);
    SaveTsv = fullfile(SaveDir, ResultTableTsv);
    save(SaveFile, 'ResultsTable');
    xASL_tsvWrite(ResultsTable, SaveTsv);
end