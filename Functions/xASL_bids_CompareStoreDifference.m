function [identical,differences,dn] = xASL_bids_CompareStoreDifference(bPrintReport,differences,dn,allFiles,iFile)
%xASL_bids_CompareStoreDifference Store the difference found in a TEXT file
%
% FORMAT: [identical,differences,dn] = xASL_bids_CompareStoreDifference(bPrintReport,differences,dn,allFiles,iFile)
%
% INPUT:
%         ...
%
% OUTPUT:
%         ...
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      Store the difference found in a TEXT file.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          ...
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


