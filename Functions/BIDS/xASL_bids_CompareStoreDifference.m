function [identical,differences,dn] = xASL_bids_CompareStoreDifference(bPrintReport,differences,dn,allFiles,iFile)
%xASL_bids_CompareStoreDifference Store the difference found in a TEXT file
%
% FORMAT: [identical,differences,dn] = xASL_bids_CompareStoreDifference(bPrintReport,differences,dn,allFiles,iFile)
%
% INPUT:
%         bPrintReport  - Print difference in file content (BOOLEAN, REQUIRED)
%         differences   - Array containing the individual differenecs (CELL ARRAY, REQUIRED)
%         dn            - Difference number for the row within differences (INTEGER, REQUIRED)
%         allFiles      - Array containing the file names (CELL ARRAY, REQUIRED)
%         iFile         - File number (INTEGER, REQUIRED)
%
% OUTPUT:
%         identical     - False if datasets are not identical
%         differences   - Array containing the individual differenecs
%         dn            - Difference number for the row within differences
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      Store the difference found in a TEXT file.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          [identical,differences,dn] = xASL_bids_CompareStoreDifference(1,differences,1,allFiles,1);
%
% REFERENCES:       ...
% __________________________________
% Copyright (c) 2015-2021 ExploreASL

    if bPrintReport
        fprintf('%s:\t\t\n',allFiles{iFile});
        fprintf('           Different file content.\n');
    end
    identical = false;
    % Save difference
    differences{dn,1} = ['Different file content: ', allFiles{iFile}, ' '];
    dn = dn+1;

end


