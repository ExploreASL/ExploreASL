function [summaryStruct] = xASL_test_bids_validation(resultDir)
    % xASL_test_bids_validation script to run over bids-validator results to consolidate errors and warnings
    %
    % FORMAT: [summaryStruct] = xASL_test_UnitTesting([bPull])
    %
    % INPUT:        resultDir       - Location of output of the bids validator (REQUIRED)
    %
    % OUTPUT:       summaryStruct   - structure containing the bids validator results
    %
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % DESCRIPTION:  Function to summarize the output of xASL_bids_validation.sh
    %
    % EXAMPLE:      [summaryStruct] = xASL_test_bids_validation('.../FlavorDatabase/2022-12-12T16:40+01:00/');
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % Copyright 2015-2022 ExploreASL
    
    %% Admin
    if nargin < 1 || isempty(resultDir)
        warning('Input must contain a folder containing output of bids-validator');
        return;
    end

    % Administration
    TimeString = datestr(now,'yyyy-mm-dd_HH_MM');
    TimeString(end-2) = 'h';
    ResultDirToday=fullfile(resultDir,TimeString);
    xASL_adm_CreateDir(ResultDirToday);

    % Generate all jsons (DEPENDENCY ON BIDS VALIDATOR) UNTESTED
    FolderList = xASL_adm_GetFileList(resultDir,'^.+$', 'List', [0 Inf],true);
    for iList=1:length(FolderList)
        iFolder = fullfile(resultDir, FolderList{iList},'rawdataReference');
        iFile   = [fullfile(ResultDirToday, FolderList{iList}), '.json'];
        input   = ['bids-validator ', iFolder , ' --no-color --json >> ', iFile];
        xASL_system(input, true);
    end

    summaryStruct = struct('errors', struct, 'warnings', struct, 'ignored', struct ,'alldata', struct());
    JsonFileList = xASL_adm_GetFileList(ResultDirToday,'^*.json$','List',[0 Inf], false);

    for iList=1:length(JsonFileList)
        % Get the full file path
        iFile = fullfile(ResultDirToday, JsonFileList{iList});
        % Filter out the name of the file without timestamp and extension.
        flavorName = string(extractBetween(JsonFileList{iList}, 1, length(JsonFileList{iList})-5));
        % Add the matfile to the summary if it's not another summary file.
        bidsOutput = spm_jsonread(iFile);

        if ~isfield(bidsOutput, 'issues')
            warning('Json file does not contain issues, skipping json')
            continue
        end
        
        summaryStruct = xASL_adm_bids_read(summaryStruct, bidsOutput.issues, flavorName);

        % Clean up any non alphanumeric characters from flavorName as to make it a struct-acceptable name
        summaryStruct.alldata.(regexprep(flavorName,'[^a-zA-Z0-9_]', "_")) = bidsOutput;

    end
end 
    
function [resultsStruct] = xASL_adm_bids_read(resultsStruct, issues, flavor)
    resultsStruct.errors   = xASL_adm_bids_readissue((resultsStruct.errors)  ,(issues.errors)  , flavor);  
    resultsStruct.warnings = xASL_adm_bids_readissue((resultsStruct.warnings),(issues.warnings), flavor);  
    resultsStruct.ignored  = xASL_adm_bids_readissue((resultsStruct.ignored) ,(issues.ignored) , flavor);  
end

function [resultsStruct] = xASL_adm_bids_readissue(resultsStruct, output, flavor)
    for iOut=1:size(output,1)
        % if the structure is empty, add initial field. 
        if ~isfield(resultsStruct, 'key')
            resultsStruct.key = {(output(iOut).key)};
            resultsStruct.flavor = {flavor};
            resultsStruct.num = {1} ;
        % if the structure doesnt have the key value for the current error, add it. 
        elseif ~contains(resultsStruct.key,(output(iOut).key))
            iKey = size(resultsStruct.key, 2) + 1;
            resultsStruct.key(iKey) = {(output(iOut).key)};
            resultsStruct.flavor(iKey) = {flavor};
            resultsStruct.num(iKey) = {1};
        % if the structure does have the key value for the current error, append it to the list. 
        else               
            iKey = find(contains(resultsStruct.key,(output(iOut).key)));
            resultsStruct.flavor(iKey) = {flavor};
            resultsStruct.num(iKey) = {(1 + cell2mat(resultsStruct.num(iKey)))};
        end
    end
end
