function [x] = xASL_imp_CheckImportSettings(x)
%xASL_imp_CheckImportSettings Basic import checks before execution
%
% FORMAT: [x] = xASL_imp_CheckImportSettings(x)
%
% INPUT:
%   x            - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Basic import checks before execution of dcm2nii
%
% 1. Check permissions for DCM2NII
% 2. Check correct DCMNII version
% 3. Check DCM Extension Filter
% 4. Check skip subjects option
% 5. Check token definitions
% 6. Manage .nii vs .nii.gz extensions
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright 2015-2023 ExploreASL


    %% 1. Check bCheckPermissions
    if x.modules.import.settings.bCheckPermissions
        dcm2niiDir = fullfile(x.opts.MyPath, 'External', 'MRIcron');
        xASL_adm_CheckPermissions(dcm2niiDir, true); % dcm2nii needs to be executable
    end
    
    %% 2. Check dcm2nii version
    if ~isfield(x.modules.import.imPar, 'dcm2nii_version') || isempty(x.modules.import.imPar.dcm2nii_version)
		x.modules.import.imPar.dcm2nii_version = '20220720';
    end
    
    %% 3. Check DCM Extension Filter
	if ~isfield(x.modules.import.imPar,'dcmExtFilter') || isempty(x.modules.import.imPar.dcmExtFilter)
        % dcmExtFilter: the last one is because some convertors save files without extension, 
		% but there would be a dot/period before a bunch of numbers
		x.modules.import.imPar.dcmExtFilter = '^(.*\.dcm|.*\.DCM|.*\.img|.*\.IMA|[^.]+|.*\.\d*)$';
	end

	%% 4. Check SkipSubjectIfExists
	if ~isfield(x.modules.import.imPar,'SkipSubjectIfExists') || isempty(x.modules.import.imPar.SkipSubjectIfExists)
		% allows to skip existing subject folders in the temp folder, when this is set to true,
		% avoiding partly re-importing/converting dcm2niiX when processing has been partly done
		x.modules.import.imPar.SkipSubjectIfExists = false;
	else
		warning('Skipping existing subjects in temp folder...');
		fprintf('If you want to overwrite, first remove the full subject folder...');
	end

	%% 5. Check if captured groups in folderHierarchy matches with the number of tokens defined in tokenOrdering
    % Any (X|Y|Z) expression is referred to as a "captured group"
    % All captured groups defined in tokenOrdering, scanAliases, visitAliases, and sessionAliases are "tokens"
    % While it is simplest if all captured groups in folderHierarchy represent tokens (for the user and bugfixing), this is not required

	nLeftBrackets = 0;
	nRightBrackets = 0;

    % Check if the number of left brackets == number of right brackets
	for iLayer = 1:numel(x.modules.import.imPar.folderHierarchy)
		leftBracket = regexp(x.modules.import.imPar.folderHierarchy{iLayer}, '(');
		rightBracket = regexp(x.modules.import.imPar.folderHierarchy{iLayer}, ')');
		if numel(leftBracket)~=numel(rightBracket)
			error(['Unequal brackets used in folderHierarchy layer ' num2str(iLayer) ': ' x.modules.import.imPar.folderHierarchy{iLayer}]);
		end

		nLeftBrackets = nLeftBrackets + numel(leftBracket);
		nRightBrackets = nRightBrackets + numel(rightBracket);
	end

    nGroupsFolderHierarchy = nLeftBrackets;
	nGroupsTokenOrdering = max(x.modules.import.imPar.tokenOrdering);
    nTokensTokenOrdering = sum(x.modules.import.imPar.tokenOrdering>0);

    nVisitAliases = size(x.modules.import.imPar.tokenVisitAliases,1);
    nSessionAliases = size(x.modules.import.imPar.tokenSessionAliases,1);
    nScanAliases = size(x.modules.import.imPar.tokenScanAliases,1);

	% Report this
	fprintf('%s\n', [xASL_num2str(nGroupsFolderHierarchy) ' captured groups () defined in folderHierarchy']);
    fprintf('%s\n', [xASL_num2str(nGroupsTokenOrdering) ' captured groups () defined in tokenOrdering']);
	fprintf('%s\n', [xASL_num2str(nTokensTokenOrdering) ' tokens () defined in tokenOrdering']);

    fprintf('%s\n', [xASL_num2str(nVisitAliases) ' visits defined in tokenVisitAliases']);
    fprintf('%s\n', [xASL_num2str(nSessionAliases) ' sessions defined in tokenSessionAliases: ']);
    fprintf('%s\n', [xASL_num2str(nScanAliases) ' scans defined in tokenSessionAliases']);
    
    if nGroupsFolderHierarchy < nTokensTokenOrdering
		error('The number of captured groups in folderHierarchy should >= the number of tokens in tokenOrdering');
    elseif nGroupsTokenOrdering < nTokensTokenOrdering
        error('The number of captured groups in tokenOrdering should >= the number of tokens in tokenOrdering');
    end

    % These warnings can be additive, so we report them separately
	if nGroupsFolderHierarchy > nTokensTokenOrdering
		warning('The number of captured groups in folderHierarchy is higher than the number of tokens in tokenOrdering. Ensure that this is correct');
    end
    if nGroupsTokenOrdering ~= nTokensTokenOrdering
		warning('Not all captured groups in folderHierarchy are used as a token in tokenOrdering. Ensure that this is correct')
    end

    % Check visits
    if nVisitAliases > 0
        visitsAre = unique(x.modules.import.imPar.tokenVisitAliases(:,2));
        nUniqueVisits = numel(visitsAre);

        if nUniqueVisits ~= nVisitAliases
            warning('Visit definitions are used more than one time. Ensure that this is correct')
            for iVisit=1:numel(visitsAre)
                fprintf('%s\n', ['Visit ' num2str(iVisit) ' = ' visitsAre{iVisit}]);
            end
        end
    end

    % Check sessions
    if nSessionAliases > 0
        sessionsAre = unique(x.modules.import.imPar.tokenSessionAliases(:,2));
        nUniqueSessions = numel(sessionsAre);
    
        if nUniqueSessions ~= nSessionAliases
            warning('Session definitions are used more than one time. Ensure that this is correct')
            for iSession=1:numel(sessionsAre)
                fprintf('%s\n', ['Session ' num2str(iSession) ' = ' sessionsAre{iSession}]);
            end
        end
    end

    fprintf('\n');

	%% 6. Manage .nii vs .nii.gz extensions
	if ~isempty(x.modules.import.imPar.folderHierarchy)
		lastHierarchy = x.modules.import.imPar.folderHierarchy{end};
		% Remove $ first
		if strcmp(lastHierarchy(end), '$')
			lastHierarchy = lastHierarchy(1:end-1);
		end
		% Fix .nii vs .nii.gz
		if length(lastHierarchy) > 4 && strcmp(lastHierarchy(end-4:end), '\.nii')
			lastHierarchy = [lastHierarchy(1:end-5) '(\.nii|\.nii\.gz)'];
		elseif length(lastHierarchy) > 8 && strcmp(lastHierarchy(end-8:end), '\.nii\.gz')
			lastHierarchy = [lastHierarchy(1:end-9) '(\.nii|\.nii\.gz)'];
		end
		% Return $
		x.modules.import.imPar.folderHierarchy{end} = [lastHierarchy '$'];
	end

end