function x = xASL_imp_ReadSourceData(x)
%xASL_imp_ReadSourceData Read source data
%
% FORMAT: x = xASL_imp_ReadSourceData(x)
%
% INPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging (REQUIRED, STRUCT)
%
% OUTPUT:
%   x                         - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%   x.modules.import.matches  - Matched files/directories
%   x.modules.import.tokens   - Tokens of files/directories
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Read source data.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright 2015-2022 ExploreASL


    %% Here we try to fix backwards compatibility
    x.modules.import.imPar = xASL_imp_TokenBackwardsCompatibility(x.modules.import.imPar);


    %% Check if directories are files are supposed to be matched
    if x.modules.import.imPar.bMatchDirectories
        strLookFor = 'Directories';
    else
        strLookFor = 'Files';
    end


    %% Check folderHierarchy
    xASL_imp_ReadSourceData_CheckFolderHierarchy(x);
    
    
    %% Start with defining the subjects, visits, sessions (i.e. BIDS runs) and scans (i.e. ScanTypes) by listing or typing
    
    % Recursively scan the directory tree using regular exspressions at each directory level. Use ()-brackets to extract tokens
    % that will be used to identify subjects, sessions and scans. In the loop below it is possible to translate the tokens
    % to more convenient strings before using them in destination paths.
    
    % Run the regexp search function
    [x.modules.import.matches, x.modules.import.tokens] = xASL_adm_FindByRegExp(...
        x.modules.import.imPar.RawRoot, x.modules.import.imPar.folderHierarchy, 'StripRoot', true, 'Match', strLookFor,'IgnoreCase',true);
    
    % Print feedback if there are no matching files
    if isempty(x.modules.import.matches)
        % This error means that there is probably something wrong with your sourceStructure.json
        fprintf(2,'Please check your sourceStructure.json file and read the import <a href="https://exploreasl.github.io/Documentation/" rel="nofollow">documentation</a>...\n');
        error('Import of sourcedata failed (no matching files)...');
	end
    % Report missing tokenOrdering field
	if ~isfield(x.modules.import.imPar,'tokenOrdering') || isempty(x.modules.import.imPar.tokenOrdering)
		error('tokenOrdering parameter not specified or empty');
	end
	
    % Copy the columns into named vectors. This construction allows for arbitrary directory hierarchies.
    
    % Make sure to select the columns that correspond to the folder ordering defined using the regular expressions above.
    
    % Define Subjects (vSubjectIDs: cell vector with extracted subject IDs (for all visits, sessions and scans)
    x.modules.import.listsIDs.vSubjectIDs = x.modules.import.tokens(:,x.modules.import.imPar.tokenOrdering(1));


end


%% ==========================================================================================================================
%% ==========================================================================================================================
%% Check the last argument of folderHierarchy
function xASL_imp_ReadSourceData_CheckFolderHierarchy(x) % PM: this should be renamed, we specifically check the last element only here!!!

    % Get last element
    lastElement = lower(x.modules.import.imPar.folderHierarchy{end});

    % Condition for file extension
    conditionFile = '\.(dcm|ima|xml|par|rec|nii|nii\.gz)';

    % Other extension
    conditionExtension = '\.';

    % Check folderHierarchy based on bMatchDirectories
    if x.modules.import.imPar.bMatchDirectories
        % Check that there is no extension in the last folderHierachy element
        % This extension should only be there if bMatchDirectories is set to false
        if ~isempty(regexpi(lastElement,conditionFile))
           warning('folderHierarchy includes a file extension but bMatchDirectories was set to true');
        end
    elseif isempty(regexpi(lastElement,conditionFile))
        % Check for extension in last folder hierachy element
           if ~isempty(regexpi(lastElement,conditionExtension))
              warning('Unknown extension in the last element of the folder hierarchy (%s)...',lastElement);
           else
              warning('No extension used in the last element of the folder hierarchy (%s)...',lastElement);
           end
    end


end