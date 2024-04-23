function xASL_adm_SortDicomToFolders(pathDICOM, nDirLayers, bUseDCMTK, bVerbose)
%xASL_adm_SortDicomToFolders Sorts DICOMs in the pathDICOM folder to directories based on the sequence name
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
% Copyright (C) 2015-2024 ExploreASL

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
	Flist   = xASL_adm_GetFileList(pathDICOM,'^.*.(?!(xlsx|ini|json))$','FPListRec',[0 Inf]);
	for iL=1:length(Flist)
		% Track progress if verbosity is on
		if bVerbose; xASL_TrackProgress(iL, length(Flist)); end

		% Read the DICOM file - all error handling is inside this function
		tDcm = xASL_io_DcmtkRead(Flist{iL}, false, bUseDCMTK);

		% Check that the header was read and contains the basic tags
		if isempty(tDcm) || ~isfield(tDcm, 'EchoTime') || ~isfield(tDcm, 'RepetitionTime') ||...
				~isfield(tDcm, 'ImageType') || isempty(tDcm.ImageType) || ~isfield(tDcm, 'ProtocolName') || isempty(tDcm.ProtocolName)
			warning(['Incomplete DICOM header: ' Flist{iL}]);
		else
			% Always add the protocol name to the directory name
			Fname = tDcm.ProtocolName;

			% Add series description if available
			if isfield(tDcm, 'SeriesDescription') && ~isempty(tDcm.SeriesDescription) && ~strcmp(tDcm.ProtocolName, tDcm.SeriesDescription)
				Fname = [Fname '_' tDcm.SeriesDescription];
			end

			if isfield(tDcm, 'SeriesNumber') && ~isempty(tDcm.SeriesNumber)
				Fname = [Fname '_' xASL_num2str(tDcm.SeriesNumber)];
			end

			% Remove special characters and create the directory name if needed
			Fname = xASL_adm_CorrectName(Fname);
			NewDir = fullfile(pathDICOM, Fname);
			xASL_adm_CreateDir(NewDir);

			[~, Pname, Pext] = fileparts(Flist{iL});
			if strcmpi(Pext, '.ima') || strcmpi(Pext, '.dcm')
				% In case the extension is IMA or DCM, we don't change the filename
				NewFile = Flist{iL};
			else
				% In other cases, we append the extension DCM
				% Note that we have to use the previous extension as this covers the care of incorrectly detecting the extension
				% when a period (.) was in the filename
				NewFile = fullfile(NewDir, [Pname Pext '.dcm']);
			end

			% The file is moved to the correct directory
			if ~strcmp(Flist{iL}, NewFile) && exist(Flist{iL}, 'file') && ~exist(NewFile, 'file')
				xASL_Move(Flist{iL}, NewFile);
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