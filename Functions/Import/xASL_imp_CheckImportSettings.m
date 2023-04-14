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
% DESCRIPTION:    Basic import checks before execution.
%
% - Check permissions for DCM2NII
% - Get correct DCMNII version
% - Define DCM Extension Filter
% - Set default for skip subjects option
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright 2015-2023 ExploreASL


    % Check bCheckPermissions
    if x.modules.import.settings.bCheckPermissions
        dcm2niiDir = fullfile(x.opts.MyPath, 'External', 'MRIcron');
        xASL_adm_CheckPermissions(dcm2niiDir, true); % dcm2nii needs to be executable
    end
    
    if ~isfield(x.modules.import.imPar, 'dcm2nii_version') || isempty(x.modules.import.imPar.dcm2nii_version)
		x.modules.import.imPar.dcm2nii_version = '20220720';
    end
    
	if ~isfield(x.modules.import.imPar,'dcmExtFilter') || isempty(x.modules.import.imPar.dcmExtFilter)
        % dcmExtFilter: the last one is because some convertors save files without extension, 
		% but there would be a dot/period before a bunch of numbers
		x.modules.import.imPar.dcmExtFilter = '^(.*\.dcm|.*\.DCM|.*\.img|.*\.IMA|[^.]+|.*\.\d*)$';
	end

	% Check SkipSubjectIfExists
	if ~isfield(x.modules.import.imPar,'SkipSubjectIfExists') || isempty(x.modules.import.imPar.SkipSubjectIfExists)
		% allows to skip existing subject folders in the temp folder, when this is set to true,
		% avoiding partly re-importing/converting dcm2niiX when processing has been partly done
		x.modules.import.imPar.SkipSubjectIfExists = false;
	else
		warning('Skipping existing subjects in temp folder...');
		fprintf('If you want to overwrite, first remove the full subject folder...');
	end

	%% Check if folderHierarchy matches with tokenOrdering
	nLeftBrackets = 0;
	nRightBrackets = 0;

	for iLayer = 1:numel(x.modules.import.imPar.folderHierarchy)
		leftBracket = regexp(x.modules.import.imPar.folderHierarchy{iLayer}, '(');
		rightBracket = regexp(x.modules.import.imPar.folderHierarchy{iLayer}, ')');
		if numel(leftBracket)~=numel(rightBracket)
			error(['Unequal brackets used in folderHierarchy layer ' num2str(iLayer) ': ' x.modules.import.imPar.folderHierarchy{iLayer}]);
		end

		nLeftBrackets = nLeftBrackets + numel(leftBracket);
		nRightBrackets = nRightBrackets + numel(rightBracket);
	end

	nTokens = max(x.modules.import.imPar.tokenOrdering);

	% Report this
	fprintf('%s\n', [num2str(nLeftBrackets) ' tokens detected in folderHierarchy']);
	fprintf('%s\n', [num2str(nTokens) ' tokens detected in tokenOrdering']);

	if nLeftBrackets<nTokens
		error('Number of tokens in folderHierarchy defined by () should be equal or higher than the length of tokenOrdering vector');
	end

	%% Manage .nii vs .nii.gz extensions
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
   