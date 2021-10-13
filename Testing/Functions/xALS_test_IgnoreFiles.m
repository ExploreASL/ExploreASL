function flavors = xALS_test_IgnoreFiles(flavors)
%xALS_test_IgnoreFiles Remove/ignore rows from flavors.comparisonTable
%
% FORMAT: flavors = xALS_test_IgnoreFiles(flavors)
%
% INPUT:
%   flavors - struct containing flavor related fields (RECOMMENDED, STRUCT)
%
% OUTPUT:
%   flavors - struct containing flavor related fields
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Remove/ignore rows from flavors.comparisonTable.
%
% EXAMPLE:      flavors = xALS_test_IgnoreFiles(flavors);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright (c) 2015-2021 ExploreASL


    %% Define the ignore list (ASL-BIDS reports, import summary files, import logs, lock files)
    ignoreList = {  'bids_report', ...
                    'import_summary', ...
                    'xASL_module_Import_', ...
                    '.status', ...
                    fullfile('ExploreASL','lock'), ...
                    fullfile('ExploreASL','log') ...
                 };
    
    %% Check flavors struct
    if ~isfield(flavors,'comparisonTable')
        warning('Missing comparisonTable...');
        return;
    end
    
    %% Iterate over ignoreList
    for iIgnore=1:numel(ignoreList)
        % Get current item
        currentItem = ignoreList{iIgnore};
        currentItem = strrep(currentItem,filesep,['\' '\\']);
        
        % Get column from comparisonTable
        allMessages = table2array(flavors.comparisonTable(:,{'message'}));
        
        % Find the rows which are supposed to be ignored
        ignoreRowsCell = regexp(allMessages,currentItem);
        ignoreRows = zeros(numel(ignoreRowsCell),1);
        for iRow=1:numel(ignoreRowsCell)
            if ~isempty(ignoreRowsCell{iRow})
                ignoreRows(iRow,1) = iRow;
            end
        end
        ignoreRows = ignoreRows(ignoreRows~=0);
        
        % Actually remove the corresponding rows
        flavors.comparisonTable(ignoreRows,:) = [];
    end
    
end

