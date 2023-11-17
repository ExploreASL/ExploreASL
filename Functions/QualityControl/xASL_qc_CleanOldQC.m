function [x] = xASL_qc_CleanOldQC(x, bRemoveCurrentSession)
%xASL_qc_RemoveOutdatedQC Removes outdated QC parameters, loads QC from JSON and cleans the current session if necessary
%
% FORMAT: [x] = xASL_qc_RemoveOutdatedQC(x[, bRemoveCurrentSession)
%
% INPUT:
%   x 	                  - structure containing fields with all information required to run this submodule (REQUIRED)
%   bRemoveCurrentSession - Removes the current session completely (OPTIONAL, DEFAULT = false)
%
% OUTPUT:
%   x        - same as input
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function removes the old and outdated QC parameters. Specifically parameters for ASL and func which are now given in a subfield by session.
%              So parameters like x.ASL.parameter are removed and we only keep x.ASL.ASL_1.parameter (and same for func and dwi).
%              It also loads further sessions from the disc and if necessary, it cleans the current session completely.
%
%     1. Load QC for all session from QC.json and save it to x-struct
%     2. Clean QC for the current session from x-struct if bRemoveCurrentSession is true
%     3. Remove outdated tags from x-struct that do not have the session sub-field for ASL, dwi, func
% EXAMPLE: x = xASL_qc_RemoveOutdatedQC(x);
% __________________________________
% Copyright (C) 2015-2023 ExploreASL	

if nargin < 2
	error('Two parameters are required');
end

if nargin < 3 || isempty(bRemoveCurrentSession)
	bRemoveCurrentSession = false;
end

% We go through all possible scantype lists
ScanTypeList = {'func', 'ASL', 'dwi'};


%%     1. Load QC for all session from QC.json and save it to x-struct
% Load QC.json and save it to x-struct
QC_Path = fullfile(x.D.ROOT, x.SUBJECT, ['QC_collection_' x.SUBJECT '.json']);
if xASL_exist(QC_Path, 'file')
	% Load the save QC parameters
	oldOutput = xASL_io_ReadJson(QC_Path);

	% Go through all possible scantypes
	for iScanType = 1:length(ScanTypeList)
		ScanType = ScanTypeList{iScanType};
		if isfield(oldOutput, ScanType)
			% Load all subfields for that scantype
			listFields = fieldnames(oldOutput.(ScanType));
			for iField = 1:length(listFields)
				% Copy only the subfields in ScanType_N
				if ~isempty(regexp(listFields{iField}, [ScanType '_\d+'], 'once'))
					% But skip the current session unless bCompleteRun is false
					if ~strcmp(listFields{iField}, x.SESSION) || ~bRemoveCurrentSession
						x.Output.(ScanType).(listFields{iField}) = oldOutput.(ScanType).(listFields{iField});
					end
				end
			end
		end
	end
end

%%     2. Clean QC for the current session from x-struct if bRemoveCurrentSession is true
% Go through all possible scantypes
for iScanType = 1:length(ScanTypeList)
	ScanType = ScanTypeList{iScanType};

	% x.Output_im.ASL should not longer be a cell array, but should instead contain subfields.s
	if isfield(x,'Output_im') && isfield(x.Output_im, ScanType) && iscell(x.Output_im.(ScanType))
		x.Output_im = rmfield(x.Output_im, ScanType);
	end

	% If doing a complete clean of the current session
	if ~isempty(regexpi(x.SESSION, [ScanType '_\d+'], 'once'))
		if bRemoveCurrentSession
			% Clear any previous QC images
			if isfield(x,'Output_im') && isfield(x.Output_im, ScanType) && isfield(x.Output_im.(ScanType), x.SESSION)
				x.Output_im.(ScanType) = rmfield(x.Output_im.(ScanType), x.SESSION);
			end

			if isfield(x, 'Output') && isfield(x.Output, ScanType) && isfield(x.Output.(ScanType), x.SESSION)
				x.Output.(ScanType) = rmfield(x.Output.(ScanType), x.SESSION);
			end
		end
	end
	%%     3. Remove outdated tags from x-struct that do not have the session sub-field for ASL, dwi, func
	% In previous versions, we have parameters for session one only directly under x.Output.ASL (or similar for other modules) - if these are detected
	% they are removed as obsolete and a warning is issued. We only keep parameters in x.Output.ASL.ASL_1 etc

	if isfield(x.Output, ScanType)
		listFields = fields(x.Output.(ScanType));

		bIssuedWarningAboutExtraFields = false;
		for iField = 1:length(listFields)
			currentField = listFields{iField};
			if isempty(regexp(currentField, [ScanType '_\d+'], 'once'))
				x.Output.(ScanType) = rmfield(x.Output.(ScanType), currentField);
				if ~bIssuedWarningAboutExtraFields
					% Obsolete fields are present, report a warning once
					warning(['QC parameters for ' ScanType ' are now provided inside field x.Output.' ScanType '.' ScanType '_n for run n (e.g. ' ScanType '_1 ' ScanType '_2 ' ScanType '_3) instead of x.Output.' ScanType ]);
					fprintf('%s\n', ['The QC output of previous ExploreASL runs of the ASL module were deleted from x.Output.' ScanType ]);
					fprintf('%s\n', '(and will be removed from the resulting QC_*.json file)');
					bIssuedWarningAboutExtraFields = true;
				end
			end
		end
	end
end
end
