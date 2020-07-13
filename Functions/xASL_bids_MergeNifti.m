function NiftiPaths = xASL_bids_MergeNifti(NiftiPaths, seqType)
% Takes a list of NII files and concatenates 3D/4D files into a 4D sequence if possible
% FORMAT: NiftiPaths = xASL_bids_MergeNifti(NiftiPaths, seqType)
% 
% INPUT:
%   NiftiPaths    - cell containing list of strings with full paths of the files (REQUIRED)
%   seqType      - Type of the file - can be 'M0' or 'ASL' (REQUIRED)
%
% OUTPUT:
% NiftiPaths      - returns either the same list of files if nothing was done or the path to the newly created file
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function takes a list of M0 or ASL4D files and concatenates them together in a longer 4D volume if possible
%              following certain patterns: works only with 3D and 4D files; all files in the list must have the same size of the
%              first three dimensions; files are generarily sorted according to the last number in the filename and outputted
%              to M0.nii or ASL4D.nii; first JSON is taken and renamed, all other JSONs and NIIs are deleted after merging;
%              M0*_parms.m or ASL*_parms.mat is renamed to M0_parms.m or ASL4D_parms.m; M0 files are checked if the field 
%              PhaseEncodingAxis is consistent through all the volumes, if not the nothing is merged; this is applied to a generic case
%              and 3 other specific Siemens scenaria are implemented:
%              i) All NII files have two volumes, then simply concatenate according to the last number.
%              ii) Two files with a single volume each are merged according to the last number in the file name.
%              iii) Multiple files with each containing a single volume are sorted to tags ASL4D_x_x_Y and controls ASL4D_Y and merged in the order
%                   of the last number in the filename (Y) alternating the tags and controls
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2020 ExploreASL, JP

%% Admin
if nargin < 1 || isempty(NiftiPaths)
	error('Requires a list of files on the input.');
end

if nargin < 2 || isempty(seqType)
	error('Files type need to be provided as M0 or ASL');
else
	if ~strcmp(seqType,'M0') && ~strcmp(seqType,'ASL')
		error('seqType must be either M0 or ASL');
	end
end


if length(NiftiPaths)>1

	switch (seqType)
		case 'M0'
			% Merges the M0 files
			pathOut = xASL_adm_MergeNiiFiles_M0Files(NiftiPaths);
			if ~isempty(pathOut)
				NiftiPaths = {pathOut};
				fprintf('Corrected dcm2niiX output for\n');
				fprintf('%s\n', pathOut);
			end
			
		case 'ASL'
			% Merges Siemens ASL file if they have the known pattern of filenames
			pathOut = xASL_adm_MergeNiiFiles_SiemensASLFiles(NiftiPaths);
			if ~isempty(pathOut)
				NiftiPaths = {pathOut};
				fprintf('Corrected dcm2niiX output for\n');
				fprintf('%s\n', pathOut);
			end
			
			% Generic merging of ASL4D files for non-Siemens, or Siemens files with an unknown pattern
			if isempty(pathOut)
				% If the previous Siemens merging didn't merge them already
				pathOut = xASL_adm_MergeNiiFiles_AllASLFiles(NiftiPaths);
								
				if ~isempty(pathOut)
					NiftiPaths = {pathOut};
					fprintf('Corrected dcm2niiX output for\n');
					fprintf('%s\n', pathOut);
				end
			end
	end
end

end

%% Generic merging of M0 files
function pathOut = xASL_adm_MergeNiiFiles_M0Files(NiftiPaths)

bCheckConsistency = 1; % So far, no error was found and files can be concatenated
pathOut = ''; % Newly assigned path of a concatenated file
listEndNumber = zeros(length(NiftiPaths),1);

% Go through all JSONs
for iFile=1:length(NiftiPaths)
	[Fpath, Ffile, ~] = xASL_fileparts(NiftiPaths{iFile});
	jsonParms = spm_jsonread(fullfile(Fpath,[Ffile,'.json']));
	
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
		if (isempty(jsonParms) || ~isfield(jsonParms,'PhaseEncodingAxis'))
			if ~isempty(strPEAxis)
				bCheckConsistency = 0;
			end
		else
			if ~strcmp(strPEAxis,jsonParms.PhaseEncodingAxis)
				bCheckConsistency = 0;
			end
		end
		
		if (isempty(jsonParms) || ~isfield(jsonParms,'PhaseEncodingDirection'))
			if ~isempty(strPEDirection)
				bCheckConsistency = 0;
			end
		else
			if ~strcmp(strPEDirection,jsonParms.PhaseEncodingDirection)
				bCheckConsistency = 0;
			end
		end
	end
end

% Check if there's no difference in AP-PA direction, if all are the same, then start merging
if bCheckConsistency
	[~,indexSortedFile] = sort(listEndNumber);
	pathOut = xASL_adm_MergeNiiFiles_Merge(NiftiPaths,indexSortedFile,'M0');
	
	if ~isempty(pathOut)
		xASL_adm_MergeNiiFiles_RenameParms(Fpath,'M0');
		xASL_adm_MergeNiiFiles_Delete(NiftiPaths);
	end
end
end

%% Merging of Siemens ASL files that have a specific pattern of filenames
function pathOut = xASL_adm_MergeNiiFiles_SiemensASLFiles(NiftiPaths)
		
bCheckConsistency = 1; % So far, no error was found and files can be concatenated
pathOut = ''; % Newly assigned path of a concatenated file
listEndNumber = zeros(length(NiftiPaths),1);

for iFile=1:length(NiftiPaths)
	[Fpath, Ffile, ~] = xASL_fileparts(NiftiPaths{iFile});
	jsonParms = spm_jsonread(fullfile(Fpath,[Ffile,'.json']));
	
	% Make sure that we deal with Siemens files
	if isempty(jsonParms) || ~isfield(jsonParms,'Manufacturer') || isempty(strfind(jsonParms.Manufacturer,'Siemens'))
		bCheckConsistency = 0;
	end
	
	% Get the size of the 4th dimension
	dimTime = xASL_io_ReadNifti(NiftiPaths{iFile});
	if length(dimTime.dat.dim) > 3
		dimTime = dimTime.dat.dim(4);
	else
		dimTime = 1;
	end
	
	% Compare the PEAxis and check that it is the same for all of the M0-files
	if iFile == 1
		checkTimeDim = dimTime;
	else
		if checkTimeDim ~= dimTime
			bCheckConsistency = 0;
		end
	end
end

% Conditions passed for merging Siemens files
if bCheckConsistency
	if checkTimeDim == 2
		% Alternating scenario C/L - merge them in the order of the last number in the filename
		for iFile=1:length(NiftiPaths)
			[Fpath, Ffile, ~] = xASL_fileparts(NiftiPaths{iFile});
			
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
			[~,indexSortedFile] = sort(listEndNumber);
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
			[~,indexSortedFileTag] = sort(listEndNumber);
			indexSortedFileTag = indexSortedFileTag(listTag(indexSortedFileTag) == 1);
			[~,indexSortedFileControl] = sort(listEndNumber);
			indexSortedFileControl = indexSortedFileControl(listTag(indexSortedFileControl) == 0);
			if length(indexSortedFileTag) == length(indexSortedFileControl)
				indexSortedFile(1:2:length(NiftiPaths)) = indexSortedFileTag;
				indexSortedFile(2:2:length(NiftiPaths)) = indexSortedFileControl;
			else
				bCheckConsistency = 0;
			end
		end
	else
		% The Siemens specific scenarios work only with single or two volumes per file
		% All others are skipped and can be merged according to the generic criteria.
		bCheckConsistency = 0;
	end
end

% Correctly identified one of the Siemens scenarios so can merge the files in the indexSortedFile order
if bCheckConsistency
	pathOut = xASL_adm_MergeNiiFiles_Merge(NiftiPaths,indexSortedFile,'ASL4D');
	
	if ~isempty(pathOut)
		xASL_adm_MergeNiiFiles_RenameParms(Fpath,'ASL4D');
		xASL_adm_MergeNiiFiles_Delete(NiftiPaths);
	end
end
end

%% Merging of any ASL files
function pathOut = xASL_adm_MergeNiiFiles_AllASLFiles(NiftiPaths)
% First rename the NIfTI and JSON files to 4 digit format & sort them
listEndNumber = zeros(length(NiftiPaths),1);
for iFile=1:length(NiftiPaths)
	[Fpath, Ffile, ~] = xASL_fileparts(NiftiPaths{iFile});
	
	% List the end number from the file name
	[iStart, iEnd] = regexp(Ffile,'\d*$');
	listEndNumber(iFile) = str2double(Ffile(iStart:iEnd));
end
[~,indexSortedFile] = sort(listEndNumber);
pathOut = xASL_adm_MergeNiiFiles_Merge(NiftiPaths,indexSortedFile,'ASL4D');

if ~isempty(pathOut)
	xASL_adm_MergeNiiFiles_RenameParms(Fpath,'ASL4D');
	xASL_adm_MergeNiiFiles_Delete(NiftiPaths);
end
end

%% Merges all the NiftiPaths in the order given by indexSortedFile and save to the pathMerged
% Save also the first JSON to pathMerged.JSON
function pathOut = xASL_adm_MergeNiiFiles_Merge(NiftiPaths,indexSortedFile,nameMerged)

bStatus = 1;
pathOut = '';
firstJSON = '';
% Start loading all the files
for iFile=1:length(NiftiPaths)
	tempIM = xASL_io_Nifti2Im(NiftiPaths{indexSortedFile(iFile)});
	
	% Merging only 3D and 4D files
	if length(size(tempIM))>4
		error('Dimensionality incorrect for this ASL NIfTI file');
	end
	
	if iFile==1
		% Get the size of the first file
		IM = tempIM;
		sizeFirst = size(IM);
		sizeFirst = sizeFirst(1:3);
	else
		% Compare size of the new file and if similar to the previous than concatenate, otherwise report an error
		sizeNew   = size(tempIM);
		sizeNew   = sizeNew(1:3);
		if isequal(sizeNew, sizeFirst)
			IM(:,:,:,end+1:end+size(tempIM,4)) = tempIM;
		else
			bStatus = 0;
		end
	end
	
	% Check for the path to JSON if existing and keep only the first existing JSON
	if isempty(firstJSON)
		[Fpath, Ffile, ~] = xASL_fileparts(NiftiPaths{indexSortedFile(iFile)});
		pathJSON = fullfile(Fpath,[Ffile '.json']);
		if exist(pathJSON,'file')
			firstJSON = pathJSON;
		end
	end
end

% If at the end and all went well
if bStatus
	fprintf('Warning: concatenating multiple NIfTIs & jsons as output from dcm2niiX\n');
	% Save the concatenated file to a given name
	pathOut = fullfile(Fpath,[nameMerged '.nii']);
	xASL_io_SaveNifti(NiftiPaths{indexSortedFile(1)}, pathOut, IM, [], 0);
	% Copy the first JSON to this name
	if ~isempty(firstJSON)
		xASL_Copy(firstJSON,fullfile(Fpath,[nameMerged '.json']),1);
	end
else
	fprintf('Warning: Cannot concatenate multiple NIfTIs & jsons as output from dcm2niiX\n');
end


end

%%
% Delete all the NiftiPaths and their associated JSONs
function xASL_adm_MergeNiiFiles_Delete(NiftiPaths)

for iFile=1:length(NiftiPaths)
	[Fpath, Ffile, ~] = xASL_fileparts(NiftiPaths{iFile});
	
	xASL_delete(NiftiPaths{iFile});
	
	pathJSON = fullfile(Fpath,[Ffile '.json']);
	% Delete JSONs
	if exist(pathJSON,'file')
		xASL_delete(pathJSON);
	end
end
end

%%
% Find file called *_parms.m in the given directory and shorten to the given name
function xASL_adm_MergeNiiFiles_RenameParms(Fpath,Fname)

FileList = xASL_adm_GetFileList(Fpath,'^*_parms.mat$','List',[],false);
if ~isempty(FileList)
	xASL_Move(fullfile(Fpath,FileList{1}),fullfile(Fpath,[Fname '_parms.mat']));
end
end