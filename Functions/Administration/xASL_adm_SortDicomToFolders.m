function xASL_adm_SortDicomToFolders(pathDICOM, nDirLayers, bUseDCMTK, bVerbose)
%xASL_adm_SortDicomToFolders Sorts DICOMs in the pathDICOM folder to directories based on the sequence name as defined in DICOM tags "SeriesDescription" or "ProtocolName". 
% Note that these are most likely the names set in the protocol on the scanner, or changed by the radiographics.
%
% FORMAT: xASL_adm_SortDicomToFolders(pathDICOM [, bUseDCMTK, bVerbose])
%
% INPUT:
%   pathDICOM    - Path to the folder with all DICOMs in directories and subdirectories (REQUIRED, STRING)
%   nDirLayers   - Number of directory layers to recursively check without renaming. For n=0, we start renaming all subdirectories at pathDICOM, 
%                  for n>1, we re-run the same function with n-1 parameter (OPTIONAL, INTEGER, DEFAULT = 0)
%   bUseDCMTK    - Use DCMTK reading, when false use SPM (OPTIONAL, BOOLEAN, DEFAULT = true)
%   bVerbose     - Verbose (OPTIONAL, BOOLEAN, DEFAULT = true)
%
% OUTPUT: It sorts all files to subdiretories and deletes all empty directories
%
% DESCRIPTION: This function sorts DICOM files into directories according to their ProtocolName/SeriesDescrption/AcquisitionNumber. 
%              It keeps the original filename, but aadds .dcm extension if needed.
%              We potentially want to run the function on a directory with several subject directories which will be intact and only their contens will be transformed.
%              The input parameter nDirLayers makes sure that the given number of directories is intact.
%
% EXAMPLE: 
%     xASL_adm_SortDicomToFolders('tmp/DICOM', [], 0)
% __________________________________
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________


if nargin<2 || isempty(nDirLayers)
	nDirLayers = 0;
end

if nargin<3 || isempty(bUseDCMTK)
    bUseDCMTK = true; % DCMTK is used by default as it is faster
end

if nargin<4 || isempty(bVerbose)
    bVerbose = true;
end

if nDirLayers > 0
	Dlist = xASL_adm_GetFsList(pathDICOM, '^.*$', 1, 0, 0, [0 Inf]);
	for iL=3:length(Dlist)
		xASL_adm_SortDicomToFolders(fullfile(pathDICOM, Dlist{iL}), nDirLayers-1, bUseDCMTK, bVerbose);
	end
else
	% Gather a list of all files in the directory
	Flist   = xASL_adm_GetFileList(pathDICOM,'.*','FPListRec',[0 Inf]);
	for iL=1:length(Flist)
		[~, ~, Fext] = xASL_fileparts(Flist{iL});
		if isempty(regexp(Fext, '(xlsx|ini|json|DS_Store)'))
			% Track progress if verbosity is on
			if bVerbose; xASL_TrackProgress(iL, length(Flist)); end

			% Read the DICOM file - all error handling is inside this function
			tDcm = xASL_io_DcmtkRead(Flist{iL}, false, bUseDCMTK, true);

			hasProtocolName = isfield(tDcm, 'ProtocolName') && ~isempty(tDcm.ProtocolName);
			hasSeriesDescription = isfield(tDcm, 'SeriesDescription') && ~isempty(tDcm.SeriesDescription);
			hasSeriesNumber = isfield(tDcm, 'SeriesNumber') && ~isempty(tDcm.SeriesNumber);
			hasSequenceName = isfield(tDcm, 'SequenceName') && ~isempty(tDcm.SequenceName);

			% Check that the header was read and contains the basic tags
			if isempty(tDcm)
                warning(['Empty DICOM header, skipping: ' Flist{iL}]);
			elseif ~hasProtocolName && ~hasSeriesDescription && ~hasSeriesNumber && ~hasSequenceName
                warning(['DICOM header without ProtocolName, SeriesDescription, SeriesNumber, or SequenceName. Skipping: ' Flist{iL}]);
			else
                % Manage directory name
                % Priority:
                % 1. ProtocolName
				% 2. SequenceName
                % 3. SeriesDescription
                % 4. SeriesNumber

                Fname = [];

                % Always add the protocol name to the directory name
                if hasProtocolName
                    Fname = tDcm.ProtocolName;
                end
                
				% Add sequence name if available
				if hasSequenceName
                    if ~isempty(Fname)
                        % if ProtocolName was available, was append SeriesDescription if it differs from ProtocolName
					    Fname = [Fname '_' tDcm.SequenceName];
                    else
                        % Only SeriesDescription is also fine
                        Fname = tDcm.SequenceName;
                    end
				end

				% Add series description if available
				if hasSeriesDescription
                    if hasProtocolName && ~strcmpi(tDcm.ProtocolName, tDcm.SeriesDescription)
                        % if ProtocolName was available, was append SeriesDescription if it differs from ProtocolName
					    Fname = [Fname '_' tDcm.SeriesDescription];
                    else
                        % Only SeriesDescription is also fine
                        Fname = tDcm.SeriesDescription;
                    end
				end                

                if isempty(Fname)
                    warning(['ProtocolName, SeriesDescription, and SequenceName missing: ' Flist{iL}]);
                end

                % Add SeriesNumber if available
				if hasSeriesNumber
					if ~isempty(Fname)
						Fname = [Fname '_'];
					end
					Fname = [Fname xASL_num2str(tDcm.SeriesNumber)];
				end
                
                %% Potential extra warnings, can disable these to reduce verbosity
                checkFields = {'EchoTime' 'RepetitionTime' 'ImageType'};
                for iField=1:length(checkFields)
                    if ~isfield(tDcm, checkFields{iField}) || isempty(tDcm.(checkFields{iField}))
                        warning([checkFields{iField} ' missing: ' Flist{iL}]);
                    end
                end

				%% Remove special characters and create the directory name if needed
				Fname = xASL_adm_CorrectName(Fname);
				NewDir = fullfile(pathDICOM, Fname);
				xASL_adm_CreateDir(NewDir);

				[~, Pname, Pext] = fileparts(Flist{iL});
				if strcmpi(Pext, '.ima') || strcmpi(Pext, '.dcm')
					% In case the extension is IMA or DCM, we don't change the filename
					NewFile = [Pname Pext];
				else
					% In other cases, we append the extension DCM
					% Note that we have to use the previous extension as this covers the care of incorrectly detecting the extension
					% when a period (.) was in the filename
					NewFile = [Pname Pext '.dcm'];
				end

                % Add the new directory
                NewFile = fullfile(NewDir, NewFile);

				% The file is moved to the correct directory
				if ~strcmp(Flist{iL}, NewFile) && exist(Flist{iL}, 'file') && ~exist(NewFile, 'file')
					xASL_Move(Flist{iL}, NewFile);
				end
			end
		end
	end

	% List all directories in the rootPath
	Dlist = xASL_adm_GetFsList(pathDICOM, '^.*$', 1, 0, 0, [0 Inf]);
	for iD=3:length(Dlist)
		% Recursively delete them including subdirectories, but skip those that are not empty (include files or non-empty directories)
		xASL_adm_SortDicomToFolders_RecursiveDirectoryDelete(pathDICOM, Dlist{iD});
	end
end
end

%% Recursively delete directories
function xASL_adm_SortDicomToFolders_RecursiveDirectoryDelete(pathROOT, dirName)
pathCurrent = fullfile(pathROOT, dirName);
% First list all subdirectories and recursively delete them
Dlist = xASL_adm_GetFsList(pathCurrent, '^.*$', 1, 0, 0, [0 Inf]);
for iD=3:length(Dlist)
	xASL_adm_SortDicomToFolders_RecursiveDirectoryDelete(pathCurrent, Dlist{iD});
end

% Then check for a list of remaining files and subdirectories
% If empty, then delete the directory
if  isempty(xASL_adm_GetFileList(pathCurrent, '^.*', 'FPListRec',[0 Inf]))
	try
		rmdir(pathCurrent);
	end
end
end