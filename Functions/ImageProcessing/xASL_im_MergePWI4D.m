function [path_PWI4D, path_PWI4D_Pop] = xASL_im_MergePWI4D(x)
% xASL_im_MergePWI4D merges several PWIs based on the list in x.modules.asl.sessionsToMerge
%
% FORMAT: [path_PWI4D, path_PWI4D_Pop] = xASL_im_MergePWI4D(x)
%
% INPUT:
%   x                   - structure containing fields with all information required to run this submodule (REQUIRED)
%
% OUTPUT: 
%   path_PWI4D          - Path to the new merged NIfTI file
%   path_PWI4D_Pop      - Path to the new merged NIfTI file in the Population directory
% OUTPUT FILES: 
%   NIfTI files (see paths above) containing the merged sessions.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function takes previously created PWI4D.nii from sessions described in the x.modules.asl.sessionsToMerge and 
%              merges them to a single NII with a predefined name. It returns this new name for both native and standard space and it also creates the 
%              appropriate JSON sidecars.
%
% EXAMPLE: 
%   [path_PWI4D, path_PWI4D_Pop] = xASL_im_MergePWI4D(x)
% __________________________________
% Copyright (C) 2015-2024 ExploreASL

% 0 Admin
if nargin<1 || isempty(x)
	error('Empty input parameters');
end

fprintf(['Concatenating PWI4D for subject ' x.SUBJECT ' for sessions: ']);
% Do this for native and standard space
pathPWI4D = {x.P.Path_PWI4D x.P.Pop_Path_PWI4D};

% Name of the merged session
sessionName = [];

% Print the name of the session being added
for iSession = 1:numel(x.modules.asl.sessionsToMerge)
	fprintf([' ' x.modules.asl.sessionsToMerge{iSession}]);
	sessionName = [sessionName '_' x.modules.asl.sessionsToMerge{iSession}];
end

for iSpace = 1:2
	for iSession = 1:numel(x.modules.asl.sessionsToMerge)
		% Here we get the path of the session, by replacing the session names 
		pathCurrentPWI4D = replace(pathPWI4D{iSpace}, x.modules.asl.sessionsToMerge{end}, x.modules.asl.sessionsToMerge{iSession});

		% Load NII and JSON
		if xASL_exist(pathCurrentPWI4D, 'file')
			[imPWI4Dcurrent, jsonPWI4Dcurrent] = xASL_io_Nifti2Im(pathCurrentPWI4D, [], [], false); % Don't even bother converting to Legacy before merging
		else
			% If one of the sessions are missing, then we issue an error
			error(['Cannot concatenate all sessions as session '  x.modules.asl.sessionsToMerge{iSession} ' is missing.']);
		end

		% Error with a missing JSON
		if isempty(jsonPWI4Dcurrent)
			error(['Cannot concatenate all sessions as session '  x.modules.asl.sessionsToMerge{iSession} ' does not have a JSON sidecar.']);
		end

		if iSession == 1
			% For the first session, we take both JSON and NII as they are
			% And apply the scaling, which is by default 1
			imPWI4DConcatenated = imPWI4Dcurrent / x.modules.asl.sessionsToMergeScaling(iSession); 

			% Note that we don't convert to Legacy and work with BIDS format directly
			jsonPWI4DConcatenated = struct();
		else	
			% For following sessions, concatenate NII and JSON
            % And apply the scaling, which is by default 1
			imPWI4Dcurrent = imPWI4Dcurrent / x.modules.asl.sessionsToMergeScaling(iSession);

			% Check dimensions
			if isequal(size(imPWI4Dcurrent, 1:3), size(imPWI4DConcatenated, 1:3))
				% Here we concatenate (cat) over the 4rd dimension (4), PWI (total concatenated image matrix) with the new current PWI
                imPWI4DConcatenated = cat(4, imPWI4DConcatenated, imPWI4Dcurrent);
			else
				error(['Cannot concatenate all sessions as session '  x.modules.asl.sessionsToMerge{iSession} ' has a different matrix size than the previous sessions.']);
			end
		end

		if isfield(jsonPWI4Dcurrent, 'EchoTime')
			if iSession == 1
				jsonPWI4DConcatenated.EchoTime = jsonPWI4Dcurrent.EchoTime;
			else
				jsonPWI4DConcatenated.EchoTime = [jsonPWI4DConcatenated.EchoTime; jsonPWI4Dcurrent.EchoTime];
			end
		else
			warning(['Missing EchoTime for session ' x.modules.asl.sessionsToMerge{iSession}]);
		end

		if isfield(jsonPWI4Dcurrent, 'PostLabelingDelay')
			if iSession == 1
				jsonPWI4DConcatenated.PostLabelingDelay = jsonPWI4Dcurrent.PostLabelingDelay;
			else
				jsonPWI4DConcatenated.PostLabelingDelay = [jsonPWI4DConcatenated.PostLabelingDelay; jsonPWI4Dcurrent.PostLabelingDelay];
			end
		else
			warning(['Missing PostLabelingDelay for session ' x.modules.asl.sessionsToMerge{iSession}]);
		end

		if isfield(jsonPWI4Dcurrent, 'LabelingDuration')
			if iSession == 1
				jsonPWI4DConcatenated.LabelingDuration = jsonPWI4Dcurrent.LabelingDuration;
			else
				jsonPWI4DConcatenated.LabelingDuration = [jsonPWI4DConcatenated.LabelingDuration; jsonPWI4Dcurrent.LabelingDuration];
			end
		else
			jsonPWI4DConcatenated.LabelingDuration = [];
			if ~strcmpi(x.Q.LabelingType, 'pasl')
				warning(['Missing LabelingDuration for session ' x.modules.asl.sessionsToMerge{iSession}]);
			end
		end

		% Take the following list of fields to save to the merged JSON
		% These fields are given per entire sequence and not per volume. So the fields are saved only once for the first session
		% and with following session, we check for consistence and report errors
		listFields2Merge = {'MagneticFieldStrength', 'PulseSequenceType', 'MRAcquisitionType', 'Manufacturer', 'ArterialSpinLabelingType'};

		% Go through all the fields to merge
		for iField2Merge = 1:length(listFields2Merge)
			% The field exists and is not empty
 			if isfield(jsonPWI4Dcurrent, listFields2Merge{iField2Merge}) && ~isempty(jsonPWI4Dcurrent.(listFields2Merge{iField2Merge}))
				if iSession == 1 
					% For the first session, just add the field
					jsonPWI4DConcatenated.(listFields2Merge{iField2Merge}) = jsonPWI4Dcurrent.(listFields2Merge{iField2Merge});
				elseif ~isfield(jsonPWI4DConcatenated, listFields2Merge{iField2Merge})
					% Session > 1, but the field was not added for session 1
					% We report a conflict
					error(['Cannot merge JSON field ' listFields2Merge{iField2Merge} ' as it was missing for the first session']);
				elseif ~isequal(jsonPWI4DConcatenated.(listFields2Merge{iField2Merge}), jsonPWI4Dcurrent.(listFields2Merge{iField2Merge}))
					% Otherwise (field added for first session and we are now adding more) check for consistency and report an error if non-consistent
					error(['Cannot merge JSON field ' listFields2Merge{iField2Merge} ' as it differs between sessions']);
				end
			elseif isfield(jsonPWI4DConcatenated, listFields2Merge{iField2Merge})
				% Field is missing, but it was added for the first session
				error(['Cannot merge JSON field ' listFields2Merge{iField2Merge} ' as it is missing for one of the sessions, but not for the first session']);
			end
		end
	end
	% Create the specific new name
	[fPath, fName] = xASL_fileparts(pathPWI4D{iSpace});
	pathOut = fullfile(fPath, [fName sessionName '.nii']);

	if iSpace == 1
		path_PWI4D = pathOut;
	else
		path_PWI4D_Pop = pathOut;
	end

	% Save NII with JSON - note that JSON is provided in BIDS so doesn't need to be converted from Legacy to BIDS
	xASL_io_SaveNifti(pathPWI4D{iSpace}, pathOut, imPWI4DConcatenated, [], 0, [], 0, jsonPWI4DConcatenated, false);

end
fprintf('\n');
end
