function [x] = xxASL_qc_RemoveOutdatedQC(x)
%xASL_qc_RemoveOutdatedQC Removes outdated QC parameters
%
% FORMAT: [x] = xASL_qc_RemoveOutdatedQC(x)
%
% INPUT:
%   x 	     - structure containing fields with all information required to run this submodule (REQUIRED)
%
% OUTPUT:
%   x        - same as input
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function removes the old and outdated QC parameters. Specifically parameters for ASL and func which are now given in a subfield by session
%
% EXAMPLE: x = xASL_qc_RemoveOutdatedQC(x);
% __________________________________
% Copyright (C) 2015-2023 ExploreASL	

% In previous versions, we have parameters for session one only directly under x.Output.ASL (or similar for other modules) - if these are detected
% they are removed as obsolete and a warning is issued. We only keep parameters in x.Output.ASL.ASL_1 etc

ScanTypeList = {'func', 'ASL'};
for iScanType = 1:length(ScanTypeList)
	ScanType = ScanTypeList{iScanType};

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
