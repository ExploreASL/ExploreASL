function [differences,identical,dn] = xASL_bids_CompareStructuresTEXT(differences,identical,bPrintReport,allFiles,iFile,dn,currentFileA,currentFileB)
%xASL_bids_CompareStructuresTEXT Compare TEXT files
%
% FORMAT: [differences,identical,dn] = xASL_bids_CompareStructuresTEXT(differences,identical,bPrintReport,allFiles,iFile,dn,currentFileA,currentFileB)
%
% INPUT:
%         ...
%
% OUTPUT:
%         ...
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      ...
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          ...
%
% REFERENCES:       ...
% __________________________________
% Copyright @ 2015-2021 ExploreASL


    % Read files if they exist
    if (exist(currentFileA,'file') && exist(currentFileB,'file')) % xASL_exist somehow didn't work here (again)
        % Compare text files content directly
        currentFileTextA = fileread(currentFileA);
        currentFileTextB = fileread(currentFileB);
        if ~strcmp(currentFileTextA,currentFileTextB)
            [identical,differences,dn] = xASL_bids_CompareStoreDifference(bPrintReport,differences,dn,allFiles,iFile);
        end
    end

end





