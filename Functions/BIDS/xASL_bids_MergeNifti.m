function [NiftiPaths, ASLContext] = xASL_bids_MergeNifti(NiftiPaths, seqType, niiTable)
%xASL_bids_MergeNifti Take a list of NIfTI files and concatenates 3D/4D files into a 4D sequence if possible
%
% FORMAT: NiftiPaths = xASL_bids_MergeNifti(NiftiPaths, seqType[, niiTable])
% 
% INPUT:
%   NiftiPaths - cell containing list of strings with full paths of the files (REQUIRED)
%   seqType    - Type of the file - can be 'M0' or 'ASL' (REQUIRED)
%   niiTable   - cell containing a table Filename, InstanceNumber, SeriesNumber, FileType, FilePath (OPTIONAL, DEFAULT = EMPTY)
%
% OUTPUT:
% NiftiPaths   - return either the same list of files if nothing was done or the path to the newly created file
% ASLContext   - ASL context, for example [deltam,m0] (CHAR ARRAY)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function takes a list of M0 or ASL4D files and concatenates them together in a longer 4D volume if possible
%              following certain patterns: works only with 3D and 4D files; all files in the list must have the same size of the
%              first three dimensions; files are generarily sorted according to the last number in the filename and outputted
%              to M0.nii or ASL4D.nii; first JSON is taken and renamed, all other JSONs and NIIs are deleted after merging;
%              M0*_parms.m or ASL*_parms.mat is renamed to M0_parms.m or ASL4D_parms.m; M0 files are checked if the field 
%              PhaseEncodingAxis is consistent through all the volumes, if not the nothing is merged; this is applied to a generic case
%              and 3 other specific Siemens scenarios are implemented:
%
%              - i) All NII files have two volumes, then simply concatenate according to the last number.
%              - ii) Two files with a single volume each are merged according to the last number in the file name.
%              - iii) Multiple files with each containing a single volume are sorted to tags ASL4D_x_x_Y and controls ASL4D_Y and merged in the order
%                   of the last number in the filename (Y) alternating the tags and controls
%
%              This function performs the following steps in subfunctions:
%
%              1. xASL_bids_MergeNifti_M0Files Generic merging of M0 files
%              2. xASL_bids_MergeNifti_GEASLFiles Merge GE ASL files and extract scan order from DICOM tags
%              3. xASL_bids_MergeNifti_SeriesNumber Merge ASL files by SeriesNumber if different
%              4. xASL_bids_MergeNifti_SiemensASLFiles Merge Siemens ASL files with specific filename pattern
%              5. xASL_bids_MergeNifti_AllASLFiles Merge any ASL files
%              6. xASL_bids_MergeNifti_Merge Merge NiftiPaths & save to pathMerged
%              7. xASL_bids_MergeNifti_Delete Delete NiftiPaths and associated JSONs
%              8. xASL_bids_MergeNifti_RenameParms Find *_parms.m files in directory and shorten to provided name
%
% EXAMPLE:     n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2021 ExploreASL


%% Admin
if nargin < 1 || isempty(NiftiPaths)
	error('Requires a list of files as input');
end

if nargin < 2 || isempty(seqType)
	error('Files type need to be provided as M0 or ASL');
else
	if ~strcmp(seqType, 'M0') && ~strcmp(seqType, 'ASL')
		error('seqType must be either M0 or ASL');
	end
end

if nargin < 3
	niiTable = {};
end

% Fallback
ASLContext = '';

if length(NiftiPaths)>1
	switch (seqType)
    case 'M0'
        % Merges the M0 files
        pathOut = xASL_bids_MergeNifti_M0Files(NiftiPaths);

    case 'ASL'
		% Run the GE merging procedure first, that returns an empty path if not all conditions are met
		[pathOut,ASLContext] = xASL_bids_MergeNifti_GEASLFiles(NiftiPaths);
		
		% Merge ASL files by SeriesNumber if different
		if isempty(pathOut) && ~isempty(niiTable)
			pathOut = xASL_bids_MergeNifti_SeriesNumber(NiftiPaths, niiTable);
		end
		
        % Merges Siemens ASL file if they have the known pattern of filenames
		if isempty(pathOut)
			pathOut = xASL_bids_MergeNifti_SiemensASLFiles(NiftiPaths);
		end

        % Generic merging of ASL4D files for non-Siemens, or Siemens files with an unknown pattern
        if isempty(pathOut)
            % If the previous Siemens merging didn't merge them already
            pathOut = xASL_bids_MergeNifti_AllASLFiles(NiftiPaths);
        end
	end
end

% If the merging worked, return the merged path
if ~isempty(pathOut)
	NiftiPaths = {pathOut};
end

end


%% ==========================================================================
%% ==========================================================================
function pathOut = xASL_bids_MergeNifti_M0Files(NiftiPaths)
%xASL_bids_MergeNifti_M0Files Generic merging of M0 files

bCheckConsistency = 1; % So far, no error was found and files can be concatenated
pathOut = ''; % Newly assigned path of a concatenated file
listEndNumber = zeros(length(NiftiPaths),1);

% Go through all JSONs
for iFile=1:length(NiftiPaths)
	[Fpath, Ffile] = xASL_fileparts(NiftiPaths{iFile});
	jsonParms = spm_jsonread(fullfile(Fpath, [Ffile, '.json']));
	
	% List the end number from the file name
	[iStart, iEnd] = regexp(Ffile,'\d*$');
	listEndNumber(iFile) = str2double(Ffile(iStart:iEnd));
	
	% Compare the PEAxis and PEDirection and check that it is the same for all of the M0-files
	if iFile == 1
		if isempty(jsonParms) || ~isfield(jsonParms,'PhaseEncodingAxis')
			strPEAxis = '';
		else
			strPEAxis = jsonParms.PhaseEncodingAxis;
		end
		
		if isempty(jsonParms) || ~isfield(jsonParms,'PhaseEncodingDirection')
			strPEDirection = '';
		else
			strPEDirection = jsonParms.PhaseEncodingDirection;
		end
	else
		if isempty(jsonParms) || ~isfield(jsonParms,'PhaseEncodingAxis')
			if ~isempty(strPEAxis)
				bCheckConsistency = 0;
			end
        elseif ~strcmp(strPEAxis,jsonParms.PhaseEncodingAxis)
				bCheckConsistency = 0;
		end
		
		if isempty(jsonParms) || ~isfield(jsonParms,'PhaseEncodingDirection')
			if ~isempty(strPEDirection)
				bCheckConsistency = 0;
			end
        elseif ~strcmp(strPEDirection,jsonParms.PhaseEncodingDirection)
				bCheckConsistency = 0;
		end
	end
end

% Check if there's no difference in AP-PA direction, if all are the same, then start merging
if bCheckConsistency
	[~, indexSortedFile] = sort(listEndNumber);
	pathOut = xASL_bids_MergeNifti_Merge(NiftiPaths, indexSortedFile, 'M0', 0);
	
	if ~isempty(pathOut)
		xASL_bids_MergeNifti_RenameParms(Fpath, 'M0');
		xASL_bids_MergeNifti_Delete(NiftiPaths);
		fprintf('Corrected dcm2niiX output for\n');
		fprintf('%s\n', pathOut);
	end
end


end

%% ==========================================================================
%% ==========================================================================
function [pathOut,ASLContext] = xASL_bids_MergeNifti_GEASLFiles(NiftiPaths)
%xASL_bids_MergeNifti_GEASLFiles merge any ASL files in alphabetical order, but also load and use the GE ASL 
% tags and save them to a correct ASL context
%
% Description: Uses the GE ImageType tag to sort out the files correctly and saves this order

% By default, the conversion did not work
pathOut = ''; 

% And ASLContext is empty
ASLContext = '';

% Goes through all files
for iFile=1:length(NiftiPaths)
	% For each file, finds the JSONs
	[jsonPath,jsonName,~] = fileparts(NiftiPaths{iFile});
	jsonPath = fullfile(jsonPath, [jsonName, '.json']);
		
	% Loads the JSON file
	if exist(jsonPath,'file')
		jsonPar = spm_jsonread(jsonPath);
	else
		fprintf('Warning: Non-existent JSON sidecar: %s\n',jsonPath) ;
		jsonPar = [];
	end
		
	% Finds the manufacturer of the file
	if ~isempty(jsonPar) && isfield(jsonPar,'Manufacturer')
		varManufacturer = jsonPar.Manufacturer;
	else
		varManufacturer = '';
	end
		
	% If GE is not identified or ImageType field doesn't exist, then exits
	if isempty(regexpi(varManufacturer,'GE')) || ~isfield(jsonPar,'ImageType')
		return;
	end
	
	% Starts looking for the correct image type
	imageType = xASL_bids_determineImageTypeGE(jsonPar);

	% If imageType is not identified for all scans, then skip this one
	if isempty(imageType)
		return;
	end
	
	% Save this to the ASL context
	if isempty(ASLContext)
		ASLContext = imageType;
	else
		ASLContext = [ASLContext,',',imageType];
	end
end

% Merges all the files together
pathOut = xASL_bids_MergeNifti_Merge(NiftiPaths,1:length(NiftiPaths),'ASL4D',0);

% If this worked
if ~isempty(pathOut)
	% And adds the ASLContext to the JSON
	[jsonPath,jsonName,~] = fileparts(pathOut);
		
	jsonPar = spm_jsonread(fullfile(jsonPath, [jsonName, '.json']));
	jsonPar.ASLContext = ASLContext;
	jsonPar = rmfield(jsonPar,'ImageType');
	spm_jsonwrite(fullfile(jsonPath, [jsonName, '.json']),jsonPar);
	
	% And deletes the old files
	xASL_bids_MergeNifti_RenameParms(jsonPath,'ASL4D');
	xASL_bids_MergeNifti_Delete(NiftiPaths);
	
	fprintf('Corrected dcm2niiX output for GE files and found the correct scan order:\n');
	fprintf('%s\n', pathOut);
	fprintf('%s\n', ASLContext);
end

end


%% ==========================================================================
%% ==========================================================================
function pathOut = xASL_bids_MergeNifti_AllASLFiles(NiftiPaths)
%xASL_bids_MergeNifti_AllASLFiles Merge any ASL files
%
% Description: First rename the NIfTI and JSON files to 4 digit format & sort them

listEndNumber = zeros(length(NiftiPaths),1);
for iFile=1:length(NiftiPaths)
	[Fpath, Ffile] = xASL_fileparts(NiftiPaths{iFile});
	
	% List the end number from the file name
	[iStart, iEnd] = regexp(Ffile,'\d*$');
	listEndNumber(iFile) = str2double(Ffile(iStart:iEnd));
end
[~, indexSortedFile] = sort(listEndNumber);
pathOut = xASL_bids_MergeNifti_Merge(NiftiPaths,indexSortedFile,'ASL4D',0);

if ~isempty(pathOut)
	xASL_bids_MergeNifti_RenameParms(Fpath,'ASL4D');
	xASL_bids_MergeNifti_Delete(NiftiPaths);
	fprintf('Corrected dcm2niiX output for following files:\n');
	fprintf('%s\n', pathOut);
end

end
