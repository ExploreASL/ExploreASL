function xASL_adm_SortDicomToFolders(pathDICOM, bUseDCMTK, bVerbose)
%xASL_adm_SortDicomToFolders Sorts DICOMs in the pathDICOM folder to directories based on the sequence name
%
% FORMAT: xASL_adm_SortDicomToFolders(pathDICOM [, bUseDCMTK, bVerbose])
%
% INPUT:
%   pathDICOM    - Path to the folder with all DICOMs in directories and subdirectories (REQUIRED, STRING)
%   bUseDCMTK    - Use DCMTK reading, when false use SPM (OPTIONAL, BOOLEAN, DEFAULT = true)
%   bVerbose     - Verbose (OPTIONAL, BOOLEAN, DEFAULT = true)
%
% OUTPUT: It sorts all files to subdiretories and deletes all empty directories
%
% DESCRIPTION: This function sorts DICOM files into directories according to their ProtocolName/SeriesDescrption/AcquisitionNumber. 
%              It keeps the original filename, but aadds .dcm extension if needed
%
% EXAMPLE: 
%     xASL_adm_SortDicomToFolders('tmp/DICOM', [], 0)
% __________________________________
% Copyright (C) 2015-2024 ExploreASL

if nargin<2 || isempty(bUseDCMTK)
    bUseDCMTK = true; % DCMTK is used by default as it is faster
end
if nargin<3 || isempty(bVerbose)
    bVerbose = true;
end

% Gather a list of all files in the directory
Flist   = xASL_adm_GetFileList(pathDICOM,'^.*.(?!(xlsx|ini|json))$','FPListRec',[0 Inf]);

for iL=1:length(Flist)
	% Track progress if verbosity is on
	if bVerbose; xASL_TrackProgress(iL, length(Flist)); end
    
	% Read the DICOM file - all error handling is inside this function
	tDcm = xASL_io_DcmtkRead(Flist{iL}, false, bUseDCMTK);

	% Check that the header was read and contains the basic tags
	if isempty(tDcm) || ~isfield(tDcm, 'EchoTime') || isempty(tDcm.EchoTime) || ~isfield(tDcm, 'RepetitionTime') || isempty(tDcm.RepetitionTime) ||...
	   ~isfield(tDcm, 'ImageType') || isempty(ItDcm.ImageType) || ~isfield(tDcm, 'ProtocolName') || isempty(tDcm.ProtocolName)
		warning(['Incomplete DICOM header: ' DicomPath]);
	else
		% Always add the protocol name to the directory name
		Fname = tDcm.ProtocolName;

		% Add series description if available
		if isfield(tDcm, 'SeriesDescription') && ~isempty(tDcm.SeriesDescription) && ~strcmp(tDcm.ProtocolName, tDcm.SeriesDescription)
			Fname = [Fname '_' tDcm.SeriesDescription];
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


Dlist = xASL_adm_GetFsList(pathDICOM,'^.*$',1,0,0,[0 Inf]);
for iD=3:length(Dlist)
	if  isempty(xASL_adm_GetFileList(fullfile(pathDICOM,Dlist{iD}),'^.*','FPListRec',[0 Inf]))
		try
			rmdir(fullfile(pathDICOM,Dlist{iD}));
		end
	end
end
end
