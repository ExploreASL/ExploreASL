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
%              dcm2nii outputs GEs deltaM and M0scan separately, even though they are scanned in a single sequence and should belong together
%              Siemens sometimes stores all controls and labels separately, which need to be reordered according their acquisition order
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
% Copyright 2015-2022 ExploreASL


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
        % 1. Merges the M0 files
        pathOut = xASL_bids_MergeNifti_M0Files(NiftiPaths);

    case 'ASL'
		% 2. Run the GE merging procedure first, that returns an empty path if not all conditions are met
		[pathOut,ASLContext] = xASL_bids_MergeNifti_GEASLFiles(NiftiPaths);
		
		% 3. Merge ASL files by SeriesNumber if different
		if isempty(pathOut) && ~isempty(niiTable)
			pathOut = xASL_bids_MergeNifti_SeriesNumber(NiftiPaths, niiTable);
		end
		
		% 4. Merges Philips ASL file if they have the known pattern of filenames
		if isempty(pathOut)
			pathOut = xASL_bids_MergeNifti_PhilipsMultiTE(NiftiPaths);
		end

        % 5. Merges Siemens ASL file if they have the known pattern of filenames
		if isempty(pathOut)
			pathOut = xASL_bids_MergeNifti_SiemensASLFiles(NiftiPaths);
		end

        % 6. Generic merging of ASL4D files for non-Siemens, or Siemens files with an unknown pattern
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






%% ===========================================================================================================
%% ===========================================================================================================
function pathOut = xASL_bids_MergeNifti_M0Files(NiftiPaths)
% 1. xASL_bids_MergeNifti_M0Files Generic merging of M0 files

pathOut = ''; % Newly assigned path of a concatenated file
listEndNumber = zeros(length(NiftiPaths),1);

% Go through all JSONs
for iFile=1:length(NiftiPaths)
	[Fpath, Ffile] = xASL_fileparts(NiftiPaths{iFile});
	jsonParms = xASL_io_ReadJson(fullfile(Fpath, [Ffile, '.json']));
	
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
				return;
			end
        elseif ~strcmp(strPEAxis,jsonParms.PhaseEncodingAxis)
				return;
		end
		
		if isempty(jsonParms) || ~isfield(jsonParms,'PhaseEncodingDirection')
			if ~isempty(strPEDirection)
				return;
			end
        elseif ~strcmp(strPEDirection,jsonParms.PhaseEncodingDirection)
				return;
		end
	end
end

% Check if there's no difference in AP-PA direction, if all are the same, then start merging
[~, indexSortedFile] = sort(listEndNumber);
pathOut = xASL_bids_MergeNifti_Merge(NiftiPaths, indexSortedFile, 'M0', 0);

if ~isempty(pathOut)
	xASL_bids_MergeNifti_RenameParms(Fpath, 'M0');
	xASL_bids_MergeNifti_Delete(NiftiPaths);
	fprintf('Corrected dcm2niiX output for\n');
	fprintf('%s\n', pathOut);
end

end






%% ===========================================================================================================
%% ===========================================================================================================

function [pathOut,ASLContext] = xASL_bids_MergeNifti_GEASLFiles(NiftiPaths)
% 2. xASL_bids_MergeNifti_GEASLFiles merge any ASL files in alphabetical order, but also load and use the GE ASL 
% tags and save them to a correct ASL context
%
% Description: Uses the GE ImageType tag to sort out the files correctly and saves this order

% By default, the conversion did not work
pathOut = ''; 

% And ASLContext is empty
ASLContext = '';

% We prepare priority scores for merging JSONs
priorityList = zeros(1,length(NiftiPaths));

% Goes through all files
for iFile=1:length(NiftiPaths)
	% For each file, finds the JSONs
	[jsonPath,jsonName,~] = xASL_fileparts(NiftiPaths{iFile});
	jsonPath = fullfile(jsonPath, [jsonName, '.json']);
		
	% Loads the JSON file
	if exist(jsonPath,'file')
		jsonPar = xASL_io_ReadJson(jsonPath);
	else
		fprintf('Warning: Non-existent JSON sidecar: %s\n', jsonPath) ;
		jsonPar = [];
	end
		
	% Finds the manufacturer of the file
	if ~isempty(jsonPar) && isfield(jsonPar, 'Manufacturer')
		varManufacturer = jsonPar.Manufacturer;
	else
		varManufacturer = '';
	end
		
	% If GE is not identified or ImageType field doesn't exist, then exits
	if isempty(regexpi(varManufacturer, 'GE')) || ~isfield(jsonPar, 'ImageType')
		return;
	end
	
	% Starts looking for the correct image type
	imageType = xASL_bids_determineImageTypeGE(jsonPar);

	% If imageType is not identified for all scans, then skip this one
	if isempty(imageType)
		return;
	end
	
	% Set priorities for merging JSONs. JSONs with low priority (low number) are merged first and then their tags are 
	% overwritten with JSONs with higher priority (higher number) 
	% Unknown file has 0 priority, then the order of increasing priority is cbf < m0scan < control/label < deltam
	switch lower(imageType)
		case 'cbf'
			priorityList(iFile) = 1;
		case 'm0scan'
			priorityList(iFile) = 2;
		case {'control', 'label'}
			priorityList(iFile) = 3;
		case 'deltam'
			priorityList(iFile) = 4;
		otherwise
			priorityList(iFile) = 0;
	end
	
	% Save this to the ASL context
	if isempty(ASLContext)
		ASLContext = imageType;
	else
		ASLContext = [ASLContext,',',imageType];
	end
end

% Merges all the files together
pathOut = xASL_bids_MergeNifti_Merge(NiftiPaths, 1:length(NiftiPaths), 'ASL4D', 0, priorityList);

% If this worked
if ~isempty(pathOut)
	% And adds the ASLContext to the JSON
	[jsonPath,jsonName,~] = xASL_fileparts(pathOut);
		
	jsonPar = xASL_io_ReadJson(fullfile(jsonPath, [jsonName, '.json']));
	jsonPar.ASLContext = ASLContext;
	jsonPar = rmfield(jsonPar,'ImageType');
	xASL_io_WriteJson(fullfile(jsonPath, [jsonName, '.json']),jsonPar);
	
	% And deletes the old files
	xASL_bids_MergeNifti_RenameParms(jsonPath,'ASL4D');
	xASL_bids_MergeNifti_Delete(NiftiPaths);
	
	fprintf('Corrected dcm2niiX output for GE files and found the correct scan order:\n');
	fprintf('%s\n', pathOut);
	fprintf('%s\n', ASLContext);
end

end






%% ===========================================================================================================
%% ===========================================================================================================
function pathOut = xASL_bids_MergeNifti_SeriesNumber(NiftiPaths, niiTable)
% 3. xASL_bids_MergeNifti_SeriesNumber Take a list of NIfTI files and
%concatenates 3D/4D files into a 4D sequence if possible according to SeriesNumber for niiTable
%
% FORMAT: pathOut = xASL_bids_MergeNifti_SeriesNumber(NiftiPaths, niiTable)
% 
% INPUT:
%   NiftiPaths - cell containing list of strings with full paths of the files (REQUIRED)
%   niiTable   - cell containing a table Filename, InstanceNumber, SeriesNumber, FileType, FilePath (REQUIRED)
% OUTPUT:
%   pathOut    - output paths
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Take a list of NIfTI files and concatenates 3D/4D files into a 4D sequence if possible according to SeriesNumber for niiTable
%
% EXAMPLE:     pathOut = xASL_bids_MergeNifti_SiemensASLFiles(NiftiPaths, niiTable);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2021 ExploreASL

% 0. Admin
% Nothing merged
pathOut = '';

% Checks if the niiTable exists
if nargin < 2 || isempty(niiTable)
	return;
end

% The SeriesNumber column exist and is numeric
if size(niiTable,2) < 5 || ~isnumeric(niiTable{1,5})
	return;
end

% All cells are filled 
if sum(cellfun('isempty',niiTable(:,5))) > 0
	return;
end

% Initialize the vector
vectorSeriesNumber = zeros(size(niiTable,1),1);
for iField = 1:size(niiTable,1)
	vectorSeriesNumber(iField,1) = niiTable{iField,5}(1);
end

if sum(isnan(vectorSeriesNumber)) > 0
	return;
end

% Sort and extract unique numbers
[~,indexSort,~] = unique(vectorSeriesNumber);

% In case all numbers are unique, then sort the files accordingly
if length(indexSort) == length(vectorSeriesNumber)
	pathOut = xASL_bids_MergeNifti_Merge(NiftiPaths,indexSort,'ASL4D',0);
end

end


%% ===========================================================================================================
%% ===========================================================================================================
function pathOut = xASL_bids_MergeNifti_PhilipsMultiTE(NiftiPaths)
% 4. xASL_bids_MergeNifti_PhilipsMultiTE Take a list of NIfTI files and
%concatenates 3D/4D files into a 4D sequence if possible (Philips Multi-TE)
%
% FORMAT: pathOut = xASL_bids_MergeNifti_PhilipsMultiTE(NiftiPaths)
% 
% INPUT:
%   NiftiPaths - cell containing list of strings with full paths of the files (REQUIRED)
%
% OUTPUT:
%   pathOut    - output paths
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Take a list of NIfTI files and concatenates 3D/4D files into a 4D sequence if possible (Philips multi-TE).
%
% EXAMPLE:     pathOut = xASL_bids_MergeNifti_PhilipsMultiTE(NiftiPaths);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2023 ExploreASL

pathOut = ''; % Newly assigned path of a concatenated file
numberTE = 1;
for iFile=1:length(NiftiPaths)
	[Fpath, Ffile] = xASL_fileparts(NiftiPaths{iFile});
	jsonParms = xASL_io_ReadJson(fullfile(Fpath,[Ffile,'.json']));

	% Make sure that we deal with Philips files
	if isempty(jsonParms) || ~isfield(jsonParms, 'Manufacturer') || isempty(strfind(jsonParms.Manufacturer, 'Philips'))
		return;
	end

	% Check that the EchoNumber field is defined for all and is increasing
	if ~isfield(jsonParms, 'EchoNumber') || jsonParms.EchoNumber ~= iFile
		return;
	end

	if jsonParms.EchoNumber > numberTE
		numberTE = jsonParms.EchoNumber;
	end
	if iFile == length(NiftiPaths)
		sizeASL = size(xASL_io_Nifti2Im(NiftiPaths{iFile}));
	end
end

% Conditions passed - we have multiple files with equal dimensions and increasing EchoTimeNumber
% So we load all of them, concatenate them and swap the PLD and TE order

% Concatenate 
pathOut = xASL_bids_MergeNifti_Merge(NiftiPaths, 1:length(NiftiPaths), 'ASL4D', 0);

if ~isempty(pathOut)
	xASL_bids_MergeNifti_RenameParms(Fpath, 'ASL4D');
	xASL_bids_MergeNifti_Delete(NiftiPaths);
	fprintf('Corrected dcm2niiX output for Philips files:\n');
	fprintf('%s\n', pathOut);

	% We need to switch the order between TE and PLD
	imConcatenated = xASL_io_Nifti2Im(pathOut);
	imReordered = zeros(size(imConcatenated));

	for iTE=1:numberTE
		imReordered(:,:,:, iTE + numberTE*(0:(sizeASL(4)-1))) = imConcatenated(:,:,:, (iTE-1)*sizeASL(4)+(1:sizeASL(4))); 
	end

	% Correctly reordered image is saved
	xASL_io_SaveNifti(pathOut, pathOut, imReordered);
end

end

%%%%%%%%%%%%%%%UPUPUPUP



%% ===========================================================================================================
%% ===========================================================================================================
function pathOut = xASL_bids_MergeNifti_SiemensASLFiles(NiftiPaths)
% 5. xASL_bids_MergeNifti_SiemensASLFiles Take a list of NIfTI files and
%concatenates 3D/4D files into a 4D sequence if possible (Siemens)
%
% FORMAT: pathOut = xASL_bids_MergeNifti_SiemensASLFiles(NiftiPaths)
% 
% INPUT:
%   NiftiPaths - cell containing list of strings with full paths of the files (REQUIRED)
%
% OUTPUT:
%   pathOut    - output paths
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Take a list of NIfTI files and concatenates 3D/4D files into a 4D sequence if possible (Siemens).
%
% EXAMPLE:     pathOut = xASL_bids_MergeNifti_SiemensASLFiles(NiftiPaths);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2023 ExploreASL


    %% xASL_bids_MergeNifti_SiemensASLFiles Merge Siemens ASL files with specific filename pattern
    bAlternatingControlLabel = 0; % == 1 when two files are given, one file is filled with controls and the other with labels
    pathOut = ''; % Newly assigned path of a concatenated file
    listEndNumber = zeros(length(NiftiPaths),1);

    for iFile=1:length(NiftiPaths)
        [Fpath, Ffile] = xASL_fileparts(NiftiPaths{iFile});
        jsonParms = xASL_io_ReadJson(fullfile(Fpath,[Ffile,'.json']));

        % Make sure that we deal with Siemens files
        if isempty(jsonParms) || ~isfield(jsonParms,'Manufacturer') || isempty(strfind(jsonParms.Manufacturer,'Siemens'))
			return;
        end

        % Get the size of the 4th dimension
        dimTime = size(xASL_io_Nifti2Im(NiftiPaths{iFile}),4);

        % Compare the PEAxis and check that it is the same for all of the M0-files
        if iFile == 1
            checkTimeDim = dimTime;
        else
            if checkTimeDim ~= dimTime
                return;
            end
        end
    end

    % Conditions passed for merging Siemens files
	if checkTimeDim == 2
		% Alternating scenario C/L - merge them in the order of the last number in the filename
		for iFile=1:length(NiftiPaths)
			[Fpath, Ffile] = xASL_fileparts(NiftiPaths{iFile});

			% List the end number from the file name
			[iStart, iEnd] = regexp(Ffile,'\d*$');
			listEndNumber(iFile) = str2double(Ffile(iStart:iEnd));
		end
		[~,indexSortedFile] = sort(listEndNumber);

	elseif checkTimeDim == 1
		% Each volume contains only a label and control
		if length(NiftiPaths) == 2
			% If there are only two of them then sort by the trailing number in the filename
			for iFile=1:length(NiftiPaths)
				[Fpath, Ffile, ~] = xASL_fileparts(NiftiPaths{iFile});

				% List the end number from the file name
				[iStart, iEnd] = regexp(Ffile,'\d*$');
				listEndNumber(iFile) = str2double(Ffile(iStart:iEnd));
			end
			[~, indexSortedFile] = sort(listEndNumber);
		else
			% If there are multiple, then break to tags ASL4D_x_x_Y.nii and controls ASL4D_Y.nii
			% Sort each by Y and alternate them
			% If there are only two of them then sort by the trailing number in the filename
			listTag = zeros(length(NiftiPaths),1);
			for iFile=1:length(NiftiPaths)
				[Fpath, Ffile, ~] = xASL_fileparts(NiftiPaths{iFile});

				% List the end number from the file name
				[iStart, iEnd] = regexp(Ffile,'\d*$');
				listEndNumber(iFile) = str2double(Ffile(iStart:iEnd));

				% Check for the x_x_Y pattern in the filename identifying the tags
				[iStart, iEnd] = regexp(Ffile,'\d*_\d*_\d*$');
				if ~isempty(iStart) && ~isempty(iEnd)
					listTag(iFile) = 1;
				else
					listTag(iFile) = 0;
				end

			end
			[~, indexSortedFileTag] = sort(listEndNumber);
			indexSortedFileTag = indexSortedFileTag(listTag(indexSortedFileTag) == 1);
			[~, indexSortedFileControl] = sort(listEndNumber);
			indexSortedFileControl = indexSortedFileControl(listTag(indexSortedFileControl) == 0);
			if length(indexSortedFileTag) == length(indexSortedFileControl)
				indexSortedFile(1:2:length(NiftiPaths)) = indexSortedFileTag;
				indexSortedFile(2:2:length(NiftiPaths)) = indexSortedFileControl;
			else
				return;
			end
		end
	elseif checkTimeDim > 2 && length(NiftiPaths) == 2
		% Exactly two files with longer time dimension is given - we count that one is control and the other label and we merge them
		% by alternating the volumes from both
		indexSortedFile = [1 2];
		bAlternatingControlLabel = 1;
	else
		% The Siemens specific scenarios work only with single or two volumes per file
		% And also with the case that exactly two files are given with more than a single time-dim
		% All others are skipped and can be merged according to the generic criteria.
		return;
	end
    

    % Correctly identified one of the Siemens scenarios so can merge the files in the indexSortedFile order
	pathOut = xASL_bids_MergeNifti_Merge(NiftiPaths,indexSortedFile,'ASL4D',bAlternatingControlLabel);

	if ~isempty(pathOut)
		xASL_bids_MergeNifti_RenameParms(Fpath,'ASL4D');
		xASL_bids_MergeNifti_Delete(NiftiPaths);
		fprintf('Corrected dcm2niiX output for Siemens files:\n');
		fprintf('%s\n', pathOut);
	end
    
end



%% ===========================================================================================================
%% ===========================================================================================================
function pathOut = xASL_bids_MergeNifti_AllASLFiles(NiftiPaths)
% 6. xASL_bids_MergeNifti_AllASLFiles Merge any ASL files
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
pathOut = xASL_bids_MergeNifti_Merge(NiftiPaths, indexSortedFile, 'ASL4D', 0);

if ~isempty(pathOut)
	xASL_bids_MergeNifti_RenameParms(Fpath,'ASL4D');
	xASL_bids_MergeNifti_Delete(NiftiPaths);
	fprintf('Corrected dcm2niiX output for following files:\n');
	fprintf('%s\n', pathOut);
end

end






%% ===========================================================================================================
%% ===========================================================================================================
function pathOut = xASL_bids_MergeNifti_Merge(NiftiPaths, indexSortedFile, nameMerged, bAlternatingControlLabel, priorityList)
% 6. xASL_bids_MergeNifti_Merge Merge NiftiPaths & save to pathOut
%
% FORMAT: pathOut = xASL_bids_MergeNifti_Merge(NiftiPaths, indexSortedFile, nameMerged, bAlternatingControlLabel[, priorityList])
% 
% INPUT:
%   NiftiPaths - Nifti paths
%   indexSortedFile - Index sorted file
%   nameMerged - Name merged
%   bAlternatingControlLabel - Alternating control label (BOOLEAN, REQUIRED)
%   priorityList - a vector of priorities indicating (higher value = higher priority) how to merge the JSON files (OPTIONAL, DEFAULT = equal priorities)
%
% OUTPUT:
%   pathOut    - output paths
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Merge NiftiPaths & save to pathOut.
%
% EXAMPLE:     ...
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2022 ExploreASL


    %% Track if all went well
    bStatus = 1; 
    pathOut = '';
	
	if nargin < 4
		error('Require 4 input parameters.');
	end
	
	if nargin < 5 || isempty(priorityList)
		priorityList = ones(1, numel(NiftiPaths));
	end
	
    %% Start loading all the files
    for iFile=1:length(NiftiPaths)
        tempIM = xASL_io_Nifti2Im(NiftiPaths{indexSortedFile(iFile)});

        % Merging only 3D and 4D files
        if length(size(tempIM))>4
            error('Dimensionality incorrect for this ASL NIfTI file');
        end

        % Compare size of the new file and if similar to the previous than concatenate, otherwise report an error
        if iFile == 1
            sizeFirst = size(tempIM);
            sizeFirst = sizeFirst(1:3);
        else
            sizeNew = size(tempIM);
            sizeNew = sizeNew(1:3);
            if ~isequal(sizeNew, sizeFirst)
                bStatus = 0;
            end
        end

        if bAlternatingControlLabel
            % Always interlace the two following files
            % For the first file, create the interleaved first volume
            if iFile == 1
                lengthFirst = size(tempIM,4);
                IM = zeros([sizeFirst,lengthFirst*2]);
                IM(:,:,:,1:2:end) = tempIM;
            elseif mod(iFile,2)
                % For odd files, create a new interleaved addition
                lengthFirst = size(tempIM,4);
                IM(:,:,:,end+1:end+2*lengthFirst) = zeros([sizeFirst,lengthFirst*2]);
                IM(:,:,:,end+2:2:end+2*lengthFirst) = tempIM;
            else
                % For even files - fill in the interleave spaces
                IM(:,:,:,end-2*lengthFirst+2:2:end) = tempIM;
            end
        else
            % Simply merge files in the order in which they come
            if iFile==1
                % Get the size of the first file
                IM = tempIM;
            else
                if bStatus
                    IM(:,:,:,end+1:end+size(tempIM,4)) = tempIM;
                end
            end
        end
        
	end

	%% Merge all the JSONs according to their priority
	
	% Create an order in which the JSONs are going to be read
	[~, priorityListIndex] = sort(priorityList);
	outputJSON = struct();
	for iJSON = 1:length(NiftiPaths)
		% Go through the JSONs in increasing priority
		[Fpath, Ffile] = xASL_fileparts(NiftiPaths{priorityListIndex(iJSON)});
		if xASL_exist(fullfile(Fpath,[Ffile '.json']), 'file')
			currentJSON = xASL_io_ReadJson(fullfile(Fpath,[Ffile '.json']));
			
			% Go through all fields
			listFieldNames = fieldnames(currentJSON);
			for iFieldName = 1:length(listFieldNames)
				% Overwrite fields as we are reading JSONs with increasing priority
				fieldsDuplicityCheck = {'GELabelingDuration','InversionTime','LabelingDuration'};
				for iField = 1:length(fieldsDuplicityCheck)
					if isfield(outputJSON, fieldsDuplicityCheck{iField}) && isfield(currentJSON, fieldsDuplicityCheck{iField})
						if ~isequal(outputJSON.(fieldsDuplicityCheck{iField}), currentJSON.(fieldsDuplicityCheck{iField}))
							warning('Difference in field %s between merged JSONs', fieldsDuplicityCheck{iField});
						end
					end
				end
				
				outputJSON.(listFieldNames{iFieldName}) = currentJSON.(listFieldNames{iFieldName});
			end
		else
			warning('While merging NIfTI files, there was a JSON file missing: %s', fullfile(Fpath,[Ffile '.json']));
		end
	end
	
	
    %% Finalize merging and save
    if bStatus
        fprintf('Warning: concatenating multiple NIfTIs & jsons as output from dcm2niiX\n');
        % Save the concatenated file to a given name
        pathOut = fullfile(Fpath,[nameMerged '.nii']);
        xASL_io_SaveNifti(NiftiPaths{indexSortedFile(1)}, pathOut, IM, [], 0);
        % Special treatment for Hadamard encoded files
        EchoTimes = cell(size(NiftiPaths,2),1);
		bHadamardFME = false;
		bMultiTE = 0;
        for iFileCheck = 1:size(NiftiPaths,2)
            % Get JSON
            [jsonPathX, jsonNameX] = xASL_fileparts(NiftiPaths{iFileCheck});
            if exist(fullfile(jsonPathX, [jsonNameX '.json']),'file')
                tmpCheckJSON = xASL_io_ReadJson(fullfile(jsonPathX, [jsonNameX '.json']));
                % Check if FME Hadamard
				% Check only once, once we find it for a single volume, then continue to assume it is FME for all volumes
				if ~bHadamardFME
					bHadamardFME = xASL_imp_CheckIfFME(tmpCheckJSON, []);
				end
				
				% Check if vendor is Siemens and contains a multi-TE setup
				
				if isfield(tmpCheckJSON,'Manufacturer') && strcmpi(tmpCheckJSON.Manufacturer,'Siemens') && isfield(tmpCheckJSON,'PhoenixProtocol')
					parameterList = xASL_bids_PhoenixProtocolReader(tmpCheckJSON.PhoenixProtocol);
					bidsPar = xASL_bids_PhoenixProtocolAnalyzer(parameterList);
					if isfield(bidsPar,'EchoTime') && length(unique(bidsPar.EchoTime)) > 1
						bMultiTE = 1;
					end
				end
				
				% For FME Hadamard case, we merge echo times
				if isfield(tmpCheckJSON,'EchoTime')
					if bHadamardFME || bMultiTE
						EchoTimes{iFileCheck,1} = tmpCheckJSON.EchoTime;
					end
				end
            end
        end
        % Add echo number array if it exists
        if sum(~cellfun(@isempty,EchoTimes))~=0
            fprintf('Merging the echo numbers of the Hadamard encoded sequence...\n');
            % Sort echo numbers
            if length(indexSortedFile)==length(EchoTimes)
                EchoTimesBackUp = EchoTimes;
                for iEchoNumber=1:length(EchoTimes)
                    EchoTimes(indexSortedFile(iEchoNumber)) = EchoTimesBackUp(iEchoNumber,1);
                end
            end
            if ~issorted(cell2mat(EchoTimes))
                fprintf('Warning: echo times do not increase, resorting will be applied...\n');
                try
                    EchoTimes = sortrows(EchoTimes,1);
                catch
                    fprintf('Sorting failed...\n');
                end
			end
			
            % Write changes to JSON
			outputJSON.EchoTime = EchoTimes;
		end
		xASL_io_WriteJson(fullfile(Fpath,[nameMerged '.json']), outputJSON);
    else
        fprintf('Warning: Cannot concatenate multiple NIfTIs & jsons as output from dcm2niiX\n');
    end

end






%% ===========================================================================================================
%% ===========================================================================================================
function xASL_bids_MergeNifti_Delete(NiftiPaths)
% 7. xASL_bids_MergeNifti_Delete Delete NiftiPaths and associated JSONs
%
% FORMAT: xASL_bids_MergeNifti_Delete(NiftiPaths);
% 
% INPUT:
%   NiftiPaths - Nifti paths
%
% OUTPUT:
%   n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Delete NiftiPaths and associated JSONs.
%
% EXAMPLE:     ...
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% Delete NiftiPaths and associated JSONs
    for iFile=1:length(NiftiPaths)
        [Fpath, Ffile] = xASL_fileparts(NiftiPaths{iFile});

        xASL_delete(NiftiPaths{iFile});

        pathJSON = fullfile(Fpath,[Ffile '.json']);
        % Delete JSONs
        xASL_delete(pathJSON); % already checks if exists before deleting
    end

end





%% ===========================================================================================================
%% ===========================================================================================================
function xASL_bids_MergeNifti_RenameParms(Fpath,Fname)
% 8. xASL_bids_MergeNifti_RenameParms Find *_parms.m files in directory and shorten to provided name
%
% FORMAT: xASL_bids_MergeNifti_RenameParms(Fpath,Fname);
% 
% INPUT:
%   Fpath - File path
%   Fname - File name
%
% OUTPUT:
%   n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Find *_parms.m files in directory and shorten to provided name.
%
% EXAMPLE:     ...
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% Find *_parms.m files in directory and shorten to provided name
    FileList = xASL_adm_GetFileList(Fpath, '^.*_parms\.mat$', 'List', [], false);

    if ~isempty(FileList)
        xASL_Move(fullfile(Fpath, FileList{1}), fullfile(Fpath,[Fname '_parms.mat']));
    end


end




