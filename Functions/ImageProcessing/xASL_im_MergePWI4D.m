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
%              We merge NIfTIs and JSONs as they are, without interpretation (the parameter interpretation, conversion BIDS->Legacy etc is done later)
%
% EXAMPLE: 
%   [path_PWI4D, path_PWI4D_Pop] = xASL_im_MergePWI4D(x)
% __________________________________
% Copyright (C) 2015-2024 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________


%% ---------------------------------------------------------------------
%% 0. Admin
if nargin<1 || isempty(x)
	error('Empty input parameters');
end

fprintf(['Concatenating PWI4D for subject ' x.SUBJECT ' for sessions: ']);
% Do this for native and standard space
pathPWI4D = {x.P.Path_PWI4D x.P.Pop_Path_PWI4D};

% Name of the merged session
mergedSessionsName = [];

% Print the name of the session being added
for iSession = 1:numel(x.modules.asl.sessionsToMerge)
	fprintf([' ' x.modules.asl.sessionsToMerge{iSession}]);
	mergedSessionsName = [mergedSessionsName '_' x.modules.asl.sessionsToMerge{iSession}];
end

for iSpace = 1:2
	for iSession = 1:numel(x.modules.asl.sessionsToMerge)
        %% ---------------------------------------------------------------------
        %% 1. Load NIfTI and JSONs
		% Here we get the path of the session, by replacing the session names 
		pathCurrentPWI4D = replace(pathPWI4D{iSpace}, x.modules.asl.sessionsToMerge{end}, x.modules.asl.sessionsToMerge{iSession});

		% Load NIfTI and JSON
		if xASL_exist(pathCurrentPWI4D, 'file')
			[imPWI4Dcurrent, jsonPWI4Dcurrent] = xASL_io_Nifti2Im(pathCurrentPWI4D, [], [], false); % We don't interpret or convert variables in this function
		else
			% If one of the sessions are missing, then we issue an error
			error(['Cannot concatenate sessions: session '  x.modules.asl.sessionsToMerge{iSession} ' is missing']);
		end

		% Error with a missing JSON
		if isempty(jsonPWI4Dcurrent)
			error(['Cannot concatenate sessions: session '  x.modules.asl.sessionsToMerge{iSession} ' does not have a JSON sidecar']);
		end

        %% ---------------------------------------------------------------------
        %% 2. Merge NIfTIs
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
				error(['Cannot concatenate sessions: session '  x.modules.asl.sessionsToMerge{iSession} ' has a different matrix size']);
			end
		end

        %% ---------------------------------------------------------------------
        %% 2. Merge JSONs

		% Take the following list of fields to save to the merged JSON. 
		% Some fields are concatenated across sessions. Other fields are given per entire sequence and not per volume. So the fields are saved only once for the first session
		% and with following session, we check for consistence and report errors

		% A list of fields
		listFields2Merge = {'MagneticFieldStrength', 'PulseSequenceType', 'MRAcquisitionType', 'Manufacturer', 'ArterialSpinLabelingType', ...% These are checked but not concatenated
			                'LabelingDuration', 'PostLabelingDelay', 'EchoTime'};% These are concatenated

		% An indication of a field is merged or concatenated
		bFields2Concatenate = [0, 0, 0, 0, 0, ...
		                       1, 1, 1];

		for iField2Merge = 1:length(listFields2Merge) % Go through all the JSON fields that might need merging
			
 			if isfield(jsonPWI4Dcurrent, listFields2Merge{iField2Merge}) && ~isempty(jsonPWI4Dcurrent.(listFields2Merge{iField2Merge}))
				% The field exists and is not empty, we therefore need to check if we can merge it across sessions
				if iSession == 1 % For the first session, just add the field
					jsonPWI4DConcatenated.(listFields2Merge{iField2Merge}) = jsonPWI4Dcurrent.(listFields2Merge{iField2Merge})(:)';
				elseif ~isfield(jsonPWI4DConcatenated, listFields2Merge{iField2Merge})
					% For later session, the field should already exist for sessions 1. If not, we report an error
					error('%s %s\n %s %s %s %s %s %s','Cannot merge sessions: For subject', x.SUBJECT, 'JSON field',  listFields2Merge{iField2Merge}, 'is missing for session', x.modules.asl.sessionsToMerge{1}, 'but not for session', x.modules.asl.sessionsToMerge{iSession});
				else
					% The field seems to be consistently present/missing/empty on both sessions. The field is either check for consistency across sessions or concatenated
					if bFields2Concatenate(iField2Merge) % For vector fields, we concatenate
						jsonPWI4DConcatenated.(listFields2Merge{iField2Merge}) = [jsonPWI4DConcatenated.(listFields2Merge{iField2Merge}), jsonPWI4Dcurrent.(listFields2Merge{iField2Merge})(:)'];
					else
						% For non-vector fields, we check consistency
						if ~isequal(jsonPWI4DConcatenated.(listFields2Merge{iField2Merge}), jsonPWI4Dcurrent.(listFields2Merge{iField2Merge})(:)')
							% We report an error if non-consistent, because this makes the ASL scans incompliant
							error('%s %s\n %s %s %s %s %s %s', 'Cannot merge sessions: For subject', x.SUBJECT, 'JSON field', listFields2Merge{iField2Merge}, 'differs between sessions', x.modules.asl.sessionsToMerge{iSession}, 'and', x.modules.asl.sessionsToMerge{1});
						end
					end
				end
			elseif isfield(jsonPWI4DConcatenated, listFields2Merge{iField2Merge}) && ~isempty(jsonPWI4DConcatenated.(listFields2Merge{iField2Merge}))
				% Field is missing for the current session, but it was added for the first session
				error('%s %s\n %s %s %s %s %s %s', 'Cannot merge sessions: For subject', x.SUBJECT, 'JSON field', listFields2Merge{iField2Merge}, 'is missing for session', x.modules.asl.sessionsToMerge{iSession}, 'but not for session', x.modules.asl.sessionsToMerge{1});
			end
		end
	end

	% Create the specific new name
	[fPath, fName] = xASL_fileparts(pathPWI4D{iSpace});
	pathOut = fullfile(fPath, [fName mergedSessionsName '.nii']);

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