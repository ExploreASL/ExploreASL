function [x, matches, tokens] = xASL_imp_ReadSourceData(x,imPar)
%xASL_imp_ReadSourceData Read source data
%
% FORMAT: [x, matches, tokens] = xASL_imp_ReadSourceData(x,imPar)
%
% INPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging (REQUIRED, STRUCT)
%   imPar    - JSON file with structure with import parameters (REQUIRED, STRUCT)
%
% OUTPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%   matches  - Matched files/directories
%   tokens   - Tokens of files/directories
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Read source data.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% Check if directories are files are supposed to be matched
    if imPar.bMatchDirectories
        strLookFor = 'Directories';
    else
        strLookFor = 'Files';
    end
    
    
    %% Start with defining the subjects, visits, sessions (i.e. BIDS runs) and scans (i.e. ScanTypes) by listing or typing
    
    % Recursively scan the directory tree using regular exspressions at each directory level. Use ()-brackets to extract tokens
    % that will be used to identify subjects, sessions and scans. In the loop below it is possible to translate the tokens
    % to more convenient strings before using them in destination paths.
    [matches, tokens] = xASL_adm_FindByRegExp(imPar.RawRoot, imPar.folderHierarchy, 'StripRoot', true, 'Match', strLookFor,'IgnoreCase',true);
    
    % Print matching files
    if isempty(matches)
        warning('No matching files, skipping');
        return;
    elseif imPar.bVerbose
        fprintf('\nMatching files (#=%g):\n',length(matches));
        for iMatch=1:size(matches,1)
            fprintf('%s\n', matches{iMatch,1});
        end
    end
    
    % Copy the columns into named vectors. This construction allows for arbitrary directory hierarchies.
    % Make sure to select the columns that correspond to the folder ordering defined using the regular expressions above.
    
    % Define Subjects (so subjects may be repeated here)
    % vSubjectIDs: cell vector with extracted subject IDs (for all visits, sessions and scans; 
    x.modules.import.listsIDs.vSubjectIDs = tokens(:,imPar.tokenOrdering(1));


end



