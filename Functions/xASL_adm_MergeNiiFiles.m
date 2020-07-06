function nii_files = xASL_adm_MergeNiiFiles(nii_files, seqType)
% Takes a list of NII files and concatenates 3D/4D files into a 4D sequence if possible
% FORMAT: nii_files = xASL_adm_MergeNiiFiles(nii_files, seqType)
% 
% INPUT:
%   nii_files    - cell containing list of strings with full paths of the files (REQUIRED)
%   seqType      - Type of the file - can be 'M0' or 'ASL' (REQUIRED)
%
% OUTPUT:
% nii_files      - returns either the same list of files if nothing was done or the path to the newly created file
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
if nargin < 1 || isempty(nii_files)
	error('Requires a list of files on the input.');
end

if nargin < 2 || isempty(seqType) 
	error('Files type need to be provided as M0 or ASL');
else
	if ~strcmp(seqType,'M0') && ~strcmp(seqType,'ASL')
		error('seqType must be either M0 or ASL');
	end
end


if length(nii_files)>1
	bCheck = 1; % So far, no error was found and files can be concatenated 
	pathOut = ''; % Newly assigned path of a concatenated file
	listEndNumber = zeros(length(nii_files),1);
	switch (seqType)
		case 'M0'
			%% Generic merging of M0 files
			
			% Go through all JSONs
			for iFile=1:length(nii_files)
				[Fpath, Ffile, ~] = xASL_fileparts(nii_files{iFile});
				jsonParms = spm_jsonread(fullfile(Fpath,[Ffile,'.json']));
				
				% List the end number from the file name
				[iStart, iEnd] = regexp(Ffile,'\d*$');
				listEndNumber(iFile) = str2double(Ffile(iStart:iEnd));
				
				% Compare the PEAxis and check that it is the same for all of the M0-files
				if iFile == 1
					if isempty(jsonParms) || ~isfield(jsonParms,'PhaseEncodingAxis')
						strPEAxis = '';
					else
						strPEAxis = jsonParms.PhaseEncodingAxis;
					end
				else
					if (isempty(jsonParms) || ~isfield(jsonParms,'PhaseEncodingAxis'))
						if ~isempty(strPEAxis)
							bCheck = 0;
						end
					else
						if ~strcmp(strPEAxis,jsonParms.PhaseEncodingAxis)
							bCheck = 0;
						end
					end
				end
			end
			
			% Check if there's no difference in AP-PA direction, if all are the same, then start merging
			if bCheck
				[~,indexSortedFile] = sort(listEndNumber);
				pathOut = xASL_adm_MergeNiiFiles_Merge(nii_files,indexSortedFile,'M0');
				
				if ~isempty(pathOut)
					xASL_adm_MergeNiiFiles_RenameParms(Fpath,'M0');
					xASL_adm_MergeNiiFiles_Delete(nii_files);
					nii_files = {pathOut};
					fprintf('Corrected dcm2niiX output for\n');
					fprintf('%s\n', pathOut);
				end
			end
		case 'ASL'
			
			%% Merging of ASL4D files for Siemens using a special pattern
			for iFile=1:length(nii_files)
				[Fpath, Ffile, ~] = xASL_fileparts(nii_files{iFile});
				jsonParms = spm_jsonread(fullfile(Fpath,[Ffile,'.json']));

				% Only works with Siemens files
				if isempty(jsonParms) || ~isfield(jsonParms,'Manufacturer') || isempty(strfind(jsonParms.Manufacturer,'Siemens'))
					bCheck = 0;
				end
				
				% Get the size of the 4th dimension
				tmpNT = xASL_io_ReadNifti(nii_files{iFile});
				if length(tmpNT.dat.dim) > 3
					tmpNT = tmpNT.dat.dim(4);
				else
					tmpNT = 1;
				end
				
				% Compare the PEAxis and check that it is the same for all of the M0-files
				if iFile == 1
					checkNT = tmpNT;
				else
					if checkNT ~= tmpNT
						bCheck = 0;
					end
				end
			end
			
			% Conditions passed for merging Siemens files
			if bCheck
				if checkNT == 2
					% Alternating scenario C/L - merge them in the order of the last number in the filename
					for iFile=1:length(nii_files)
						[Fpath, Ffile, ~] = xASL_fileparts(nii_files{iFile});
						
						% List the end number from the file name
						[iStart, iEnd] = regexp(Ffile,'\d*$');
						listEndNumber(iFile) = str2double(Ffile(iStart:iEnd));
					end
					[~,indexSortedFile] = sort(listEndNumber);

				elseif checkNT == 1
					% Each volume contains only a label and control
					if length(nii_files) == 2
						% If there are only two of them then sort by the trailing number in the filename
						for iFile=1:length(nii_files)
							[Fpath, Ffile, ~] = xASL_fileparts(nii_files{iFile});
						
							% List the end number from the file name
							[iStart, iEnd] = regexp(Ffile,'\d*$');
							listEndNumber(iFile) = str2double(Ffile(iStart:iEnd));
						end
						[~,indexSortedFile] = sort(listEndNumber);
					else
						% If there are multiple, then break to tags ASL4D_x_x_Y.nii and controls ASL4D_Y.nii
						% Sort each by Y and alternate them
						% If there are only two of them then sort by the trailing number in the filename
						listTag = zeros(length(nii_files),1);
						for iFile=1:length(nii_files)
							[Fpath, Ffile, ~] = xASL_fileparts(nii_files{iFile});
						
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
							indexSortedFile(1:2:length(nii_files)) = indexSortedFileTag;
							indexSortedFile(2:2:length(nii_files)) = indexSortedFileControl;
						else
							bCheck = 0;
						end
					end
				else
					% The Siemens specific scenarios work only with single or two volumes per file
					% All others are skipped and can be merged according to the generic criteria.
					bCheck = 0;
				end
			end
				
			% Correctly identified one of the Siemens scenarios so can merge the files in the indexSortedFile order
			if bCheck
				pathOut = xASL_adm_MergeNiiFiles_Merge(nii_files,indexSortedFile,'ASL4D');
				
				if ~isempty(pathOut)
					xASL_adm_MergeNiiFiles_RenameParms(Fpath,'ASL4D');
					xASL_adm_MergeNiiFiles_Delete(nii_files);
					nii_files = {pathOut};
					fprintf('Corrected dcm2niiX output for\n');
					fprintf('%s\n', pathOut);
				end
			end
			
			%% Generic merging of ASL4D files for non-Siemens or failed-Siemens
			if isempty(pathOut)
				% First rename the NIfTI and JSON files to 4 digit format & sort them
				listEndNumber = zeros(length(nii_files),1);
				for iFile=1:length(nii_files)
					[Fpath, Ffile, ~] = xASL_fileparts(nii_files{iFile});
					
					% List the end number from the file name
					[iStart, iEnd] = regexp(Ffile,'\d*$');
					listEndNumber(iFile) = str2double(Ffile(iStart:iEnd));
				end
				[~,indexSortedFile] = sort(listEndNumber);
				pathOut = xASL_adm_MergeNiiFiles_Merge(nii_files,indexSortedFile,'ASL4D');
				
				if ~isempty(pathOut)
					xASL_adm_MergeNiiFiles_RenameParms(Fpath,'ASL4D');
					xASL_adm_MergeNiiFiles_Delete(nii_files);
					nii_files = {pathOut};
					fprintf('Corrected dcm2niiX output for\n');
					fprintf('%s\n', pathOut);
				end
			end
	end
end

end

%%
% Merges all the nii_files in the order given by indexSortedFile and save to the pathMerged
% Save also the first JSON to pathMerged.JSON
function pathOut = xASL_adm_MergeNiiFiles_Merge(nii_files,indexSortedFile,nameMerged)

bStatus = 1;
pathOut = '';
firstJSON = '';
% Start loading all the files
for iFile=1:length(nii_files)
	tempIM = xASL_io_Nifti2Im(nii_files{indexSortedFile(iFile)});
	
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
		[Fpath, Ffile, ~] = xASL_fileparts(nii_files{indexSortedFile(iFile)});
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
	xASL_io_SaveNifti(nii_files{indexSortedFile(1)}, pathOut, IM, [], 0);
	% Copy the first JSON to this name
	if ~isempty(firstJSON)
		xASL_Copy(firstJSON,fullfile(Fpath,[nameMerged '.json']),1);
	end
else
	fprintf('Warning: Cannot concatenate multiple NIfTIs & jsons as output from dcm2niiX\n');
end
	

end

%%
% Delete all the nii_files and their associated JSONs
function xASL_adm_MergeNiiFiles_Delete(nii_files)

for iFile=1:length(nii_files)
	[Fpath, Ffile, ~] = xASL_fileparts(nii_files{iFile});
	
	xASL_delete(nii_files{iFile});
	
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

fList = xASL_adm_GetFileList(Fpath,'^*_parms.m$','List',[],false);
if ~isempty(fList)
	xASL_Move(fullfile(Fpath,fList{1}),fullfile(Fpath,[Fname '_parms.m']));
end
end