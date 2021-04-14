function xASL_bids_DCM2NII(imPar, bCopySingleDicoms, bUseDCMTK, bCheckPermissions, bClone2Source,x)
%xASL_bids_DCM2NII Run the dcm2nii part of the import.
%
% FORMAT: xASL_bids_DCM2NII(imPar, bCopySingleDicoms, bUseDCMTK, bCheckPermissions, bClone2Source,x)
% 
% INPUT:
%   imPar              - JSON file with structure with import parameters (REQUIRED, STRUCT)
%   bCopySingleDicoms  - Copy Single Dicoms (REQUIRED, BOOLEAN)
%   bUseDCMTK          - Use DCMTK for DCM2NII import (REQUIRED, BOOLEAN)
%   bCheckPermissions  - Check user access rights (REQUIRED, BOOLEAN)
%   bClone2Source      - Clone to source (REQUIRED, BOOLEAN)
%   x                  - ExploreASL x structure (REQUIRED, STRUCT)
%
% OUTPUT:
%   n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Run the dcm2nii part of the import.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     xASL_bids_DCM2NII(imPar, bCopySingleDicoms, bUseDCMTK, bCheckPermissions, bClone2Source,x);
% __________________________________
% Copyright 2015-2021 ExploreASL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize defaults of dcm2nii
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dcm2niiCatchedErrors = struct; % initialization
if bCheckPermissions
	dcm2niiDir = fullfile(x.MyPath, 'External', 'MRIcron');
	xASL_adm_CheckPermissions(dcm2niiDir, true); % dcm2nii needs to be executable
end
if ~isfield(imPar,'dcm2nii_version') || isempty(imPar.dcm2nii_version)
    imPar.dcm2nii_version = '20190902'; % OR for PARREC imPar.dcm2nii_version = '20101105'; THIS IS AUTOMATED BELOW
end
if ~isfield(imPar,'dcmExtFilter') || isempty(imPar.dcmExtFilter)
    imPar.dcmExtFilter = '^(.*\.dcm|.*\.img|.*\.IMA|[^.]+|.*\.\d*)$'; % the last one is because some convertors save files without extension, but there would be a dot/period before a bunch of numbers
end
	
	if ~isfield(imPar,'SkipSubjectIfExists') || isempty(imPar.SkipSubjectIfExists)
		% allows to skip existing subject folders in the analysis folder, when this is set to true,
		% avoiding partly re-importing/converting dcm2niiX when processing has been partly done
		imPar.SkipSubjectIfExists = false;
	else
		warning('Skipping existing subjects in analysis folder');
		fprintf('If you want to overwrite, first remove the full subject folder');
	end
	
	% Create the basic folder structure for sourcedata & derivative data
	if ~exist(imPar.RawRoot, 'dir')
		warning(['Couldnt find ' imPar.RawRoot ', trying to find a different folder instead...']);
		
		% find any folder except for analysis, source, raw, derivatives
		% xASL_adm_GetFileList uses regular expressions, to create a nice list of foldernames,
		% with/without FullPath (FPList), with/without recursive (FPListRec)
		% very powerful once you know how these work
		FolderNames = xASL_adm_GetFileList(fullfile(imPar.RawRoot, imPar.studyID), '^(?!(analysis|derivatives|source|sourcedata)).*$', 'FPList', [0 Inf], true);
		
		if length(FolderNames)==1
			imPar.RawRoot = FolderNames{1};
			fprintf('%s\n', ['Found ' imPar.RawRoot ' as sourcedata folder']);
		else
			error('Couldnt find a sourcedata folder, please rename one, or move other folders');
		end
	end
	
	xASL_adm_CreateDir(imPar.AnalysisRoot);
	
	if bCheckPermissions
		xASL_adm_CheckPermissions(imPar.RawRoot, false); % don"t need execution permisions
		xASL_adm_CheckPermissions(imPar.AnalysisRoot, false);  % don"t need execution permisions
	end
	
	% here we try to fix backwards compatibility, but this may break
	if length(imPar.tokenOrdering)==3 % backwards compatibility Visits
		if size(imPar.tokenOrdering,1) > 1
			% Vertical vector
			imPar.tokenOrdering = imPar.tokenOrdering';
		end
		imPar.tokenOrdering = [imPar.tokenOrdering(1) 0 imPar.tokenOrdering(2:3)]; % insert Visits(none)
	elseif length(imPar.tokenOrdering)==2 % backwards compatibility Visits & sessions
		imPar.tokenOrdering = [imPar.tokenOrdering(1) 0 0 imPar.tokenOrdering(2)]; % insert sessions & visits(none)
	end
	
	% Path to the dictionary to initialize - we need to keep track if the dictionary has been set, because
	% Dicominfo can be used despite bUSEDCMTK==1 when DCMTK fails
	pathDcmDict = fullfile(x.MyPath,'External','xASL_DICOMLibrary.txt');
	if ~bUseDCMTK
		% -----------------------------------------------------------------------------
		% Initialize dicom dictionary by appending private philips stuff to a temporary copy
		% -----------------------------------------------------------------------------
		dicomdict('set', pathDcmDict);
	end
	
	fid_summary = -1; % initialize to be able to catch errors and close if valid
	
	% change dcmnii_version for PARREC if needed
	if ~isempty(strfind(char(imPar.folderHierarchy(end)),'PAR'))
		imPar.dcm2nii_version = '20101105';
	end
	
	% redirect output to a log file
	diary_filepath = fullfile(imPar.AnalysisRoot, ['import_log_' imPar.studyID '_' datestr(now,'yyyymmdd_HHMMSS') '.txt']);
	diary(diary_filepath);
	
	if imPar.bMatchDirectories
		strLookFor = 'Directories';
	else
		strLookFor = 'Files';
	end
	
	%
	% -----------------------------------------------------------------------------
	% Start with defining the subjects, visits, sessions (i.e. BIDS runs) and scans (i.e. ScanTypes) by listing or typing
	% -----------------------------------------------------------------------------
	
	% Recursively scan the directory tree using regular exspressions at each directory level. Use ()-brackets to extract tokens
	% that will be used to identify subjects, sessions and scans. In the loop below it is possible to translate the tokens
	% to more convenient strings before using them in destination paths.
	[matches, tokens] = xASL_adm_FindByRegExp(imPar.RawRoot, imPar.folderHierarchy, 'StripRoot', true, 'Match', strLookFor,'IgnoreCase',true);
	if isempty(matches)
		warning('No matching files, skipping');
		return;
	elseif imPar.bVerbose
		fprintf('\nMatching files:\n');
		disp(matches);
		fprintf('#=%g\n',length(matches));
	end
	
	% Copy the columns into named vectors. This construction allows for arbitrary directory hierarchies.
	% Make sure to select the columns that correspond to the folder ordering defined using the regular expressions above.
	% Define Subjects
	
	% SUBJECTS
	vSubjectIDs = tokens(:,imPar.tokenOrdering(1)); % cell vector with extracted subject IDs (for all sessions and scans)
	
	% VISITS
	if imPar.tokenOrdering(2)==0
		% a zero means: no visits applicable
		bUseVisits = false;
		vVisitIDs = cellfun(@(y) '1', vSubjectIDs, 'UniformOutput', false); % each subject has a single visit
		imPar.tokenVisitAliases = {'^1$', '_1'};
	else
		bUseVisits = true;
		vVisitIDs = tokens(:,imPar.tokenOrdering(2)); % cell vector with extracted session IDs (for all subjects and scans)
	end
	
	% SESSIONS
	if imPar.tokenOrdering(3)==0
		% a zero means: no sessions applicable
		bUseSessions = false;
		vSessionIDs = cellfun(@(y) '1', vSubjectIDs, 'UniformOutput', false); % each subject-visit has a single session
		imPar.tokenSessionAliases = {'^1$', 'ASL_1'};
	else
		bUseSessions = true;
		vSessionIDs = tokens(:,imPar.tokenOrdering(3)); % cell vector with extracted session IDs (for all subjects and scans)
	end
	
	% SCANTYPES
	vScanIDs = tokens(:,imPar.tokenOrdering(4)); % cell vector with extracted scan IDs (for all subjects and sessions)
	
	% Convert the vectors to unique & sort sets by the output aliases
	subjectIDs  = sort(unique(vSubjectIDs));
	nSubjects = length(subjectIDs);
	
	visitIDs  = unique(vVisitIDs);
	% sort by output
	if length(visitIDs)>1
		for iV=1:length(visitIDs)
			IDrow(iV) = find(cellfun(@(y) strcmp(y,visitIDs{iV}), imPar.tokenVisitAliases(:,1)));
		end
		visitIDs = visitIDs(IDrow);
	end
	nVisits = length(visitIDs);
	
	sessionIDs  = sort(unique(vSessionIDs));
	nSessions = length(sessionIDs);
	
	scanIDs = sort(unique(lower(vScanIDs)));
	nScans = length(scanIDs);
	
	% VISIT NAMES
	if isempty(imPar.visitNames)
		if isempty(visitIDs)
			imPar.visitNames = cell(nVisits,1);
			for kk=1:nVisits
				imPar.visitNames{kk}=sprintf('ASL_%g',kk);
			end
		else
			imPar.visitNames = visitIDs;
		end
	end
	
	% SESSION NAMES
	% optionaly we can have human readble session names; by default they are the same as the original tokens in the path
	if isempty(imPar.sessionNames)
		if isempty(sessionIDs)
			imPar.sessionNames = cell(nSessions,1);
			for kk=1:nSessions
				imPar.sessionNames{kk}=sprintf('ASL_%g',kk);
			end
		else
			imPar.sessionNames = sessionIDs;
		end
	end
	
	% SCAN NAMES
	scanNames = scanIDs;
	
	% sanity check for missing elements
	if nSubjects==0
		error('No subjects')
	end
	if nVisits==0
		error('No visits')
	end
	if nSessions==0
		error('No sessions')
	end
	if nScans==0
		error('No scans')
	end
	
	% preallocate space for (global) counts
	converted_scans = zeros(nSubjects,nVisits,nSessions,nScans,'uint8'); % keep a count of all individual scans
	skipped_scans = zeros(nSubjects,nVisits,nSessions,nScans,'uint8'); % keep a count of all individual scans
	missing_scans = zeros(nSubjects,nVisits,nSessions,nScans,'uint8'); % keep a count of all individual scans
	
	% define a cell array for storing info for parameter summary file
	xASL_adm_CreateDir(imPar.AnalysisRoot);
	summary_lines = cell(nSubjects,nVisits,nSessions,nScans);
	
	% -----------------------------------------------------------------------------
	% import subject by subject, visit by visit, session by session, scan by scan
	% -----------------------------------------------------------------------------
	fprintf('%s\n', 'Running import (i.e. dcm2niiX)');
	separatorline = repmat(char('+'),1,80);
	for iSubject=1:nSubjects
		subjectID = subjectIDs{iSubject};
		
		for iVisit=1:nVisits
			visitID = visitIDs{iVisit};
			
			% convert visit ID to a suitable name
			if size(imPar.tokenVisitAliases,2)==2
				iAlias = find(~cellfun(@isempty,regexp(visitIDs{iVisit},imPar.tokenVisitAliases(:,1),'once')));
				if ~isempty(iAlias)
					imPar.visitNames{iVisit} = imPar.tokenVisitAliases{iAlias,2};
				end
			end
			
			if bUseVisits % only pad VisitID _1 _2 _3 etc if there are visits specified
				% Multiple visits is defined by the tokenVisitAliases.
				% If this is non-existing, it is set to 1, and if it does exist,
				% it will put the _1 _2 _3 etc in the folder
				% this fix allows to import a single visit from a range of
				% specified visits
				SubjDir = fullfile(imPar.AnalysisRoot, [subjectID imPar.visitNames{iVisit}]);
				% if strcmp(imPar.visitNames{iVisit},'_1') % only pad the visitID _1 _2 _3 etc if there are multiple visits
			else
				SubjDir = fullfile(imPar.AnalysisRoot, subjectID);
			end
			
			if imPar.SkipSubjectIfExists && exist(SubjDir, 'dir')
				continue; % we found the subject dir (i.e. SubjectVisit), so we skip it
				% this is ignored when imPar.SkipSubjectIfExists is set to
				% false (default)
			end
			
			fprintf('%s\nImporting subject=%s:   \n',separatorline,[subjectID imPar.visitNames{iVisit}]); % display subject-visit ID
			
			% loop through all sessions
			for iSession=1:nSessions
				sessionID = sessionIDs{iSession};
				
				% convert session ID to a suitable name
				if size(imPar.tokenSessionAliases,2)==2
					iAlias = find(~cellfun(@isempty,regexp(sessionID,imPar.tokenSessionAliases(:,1),'once')));
					if ~isempty(iAlias)
						imPar.sessionNames{iSession} = imPar.tokenSessionAliases{iAlias,2};
					end
				end
				
				for iScan=1:nScans
					scanID = scanIDs{iScan};
					summary_line = [];
					first_match = [];
					summary_lines{iSubject,iVisit,iSession,iScan} = 'n/a';
					
					if ~imPar.bVerbose % if not verbose, track % progress
						CounterN = (iSession-1)*nScans+iScan;
						CounterT = nSessions*nScans;
						xASL_TrackProgress(CounterN, CounterT);
					end
					
					% convert scan ID to a suitable name and set scan-specific parameters
					if size(imPar.tokenScanAliases,2)==2
						iAlias = find(~cellfun(@isempty,regexpi(scanID,imPar.tokenScanAliases(:,1),'once')));
						if ~isempty(iAlias)
							scanNames{iScan} = imPar.tokenScanAliases{iAlias,2};
						else
							% keep the original name
							WarningMessage = ['ExploreASL_Import: Unknown scan ID ' scanID ' found, don"t know what this is'];
							dcm2niiCatchedErrors = xASL_bids_CatchErrors('isempty(iAlias)', WarningMessage, dbstack, mfilename, pwd, scan_name, scanpath, destdir, dcm2niiCatchedErrors, imPar);
						end
					end
					scan_name = scanNames{iScan};
					
					% minimalistic feedback of where we are
					if imPar.bVerbose; fprintf('>>> Subject=%s, visit=%s, session=%s, scan=%s\n',subjectID, visitID, num2str(iSession), scan_name); end
					
					bOneScanIsEnough = false; % default
					bPutInSessionFolder = true; % by default put in session folder
					switch scan_name
						case {'ASL4D', 'M0', 'ASL4D_RevPE', 'func_bold'}
							bPutInSessionFolder = true;
						case {'T1', 'WMH_SEGM', 'FLAIR', 'T2', 'T1c'}
							bPutInSessionFolder = false;
					end
					
					if ~isempty(strfind(char(imPar.folderHierarchy(end)),'PAR'))
						imPar.dcm2nii_version = '20101105';
					end
					
					% now pick the matching one from the folder list
					iMatch = find(strcmp(vSubjectIDs,subjectID) & strcmp(vVisitIDs, xASL_adm_CorrectName(visitID,2,'_')) & strcmp(vSessionIDs,sessionID) & strcmpi(vScanIDs,scanID) ); % only get the matching session
					if isempty(iMatch)
						% only report as missing if we need a scan for each session (i.e. ASL)
						if sum(converted_scans(iSubject,iVisit,:,iScan))==0
							WarningMessage = ['Missing scan: ' [subjectID imPar.visitNames{iVisit}] ', ' num2str(iSession) ', ' scan_name];
							if imPar.bVerbose; warning(WarningMessage); end
							missing_scans(iSubject, iVisit, iSession, iScan) = 1;
						end
						
						summary_lines{iSubject, iVisit, iSession, iScan} = summary_line;
						continue;
						warning('Dont forget to comment continue here for debugging');
					end
					
					% determine input and output paths
					bSkipThisOne = false;
					branch = matches{iMatch};
					scanpath = fullfile(imPar.RawRoot,branch);
					
					if ~isempty(strfind(scanNames{iScan}, 'ASL4D')) || ~isempty(strfind(scanNames{iScan}, 'M0'))
						session_name = ['ASL_' num2str(iSession)];
					elseif ~isempty(strfind(scanNames{iScan}, 'DSC4D'))
						session_name = ['DSC_' num2str(iSession)];
					else
						session_name = [scanNames{iScan} '_' num2str(iSession)]; % Allow multiple ScanTypes for sessions
					end
					
					if bPutInSessionFolder
						destdir = fullfile(SubjDir, session_name);
					else % put in subject folder instead of session folder
						destdir = SubjDir;
					end
					
					if bOneScanIsEnough && sum(converted_scans(iSubject,iVisit,:,iScan))~=0
						% one scan is enough, so skip this one if there was already a scan converted of this type (i.e. T1)
						if imPar.bVerbose; fprintf('Skipping scan: %s, %s, %s\n',[subjectID imPar.visitNames{iVisit}],session_name,scan_name); end
						bSkipThisOne = true;
						destdir = []; % just in case
					end
					
					% start the conversion if this scan should not be skipped
					if bSkipThisOne
						summary_line = sprintf(',"skipped",,,,,,,,');
						skipped_scans(iSubject, iVisit, iSession, iScan) = 1;
					else
						nii_files = {};
						xASL_adm_CreateDir(destdir);
						
						% check if we have a nii(gz) file, or something that needs to be converted (parrec/dicom)
						if ~exist(scanpath, 'dir') && ~isempty(regexpi(scanpath,'(\.nii|\.nii\.gz)$'))
							% we found a NIfTI file
							% check if output exists
							first_match = fullfile(destdir, [scan_name '.nii']);
							if imPar.bOverwrite || ~xASL_exist(first_match,'file')
								[~, fname, fext] = fileparts(scanpath);
								destfile = fullfile(destdir, [fname fext]); % will be renamed later
								xASL_Copy(scanpath, destfile, imPar.bOverwrite, imPar.bVerbose);
								% gunzip if required
								destfile = xASL_adm_UnzipNifti(destfile);
								xASL_Move(destfile, first_match, imPar.bOverwrite, imPar.bVerbose);
							end
							nii_files{1} = first_match;
						else % we found dicom files
							% -----------------------------------------------------------------------------
							% start the conversion. Note that the dicom filter is only in effect when a directory is specified as input.
							% -----------------------------------------------------------------------------
							try
								[nii_files, scan_name, first_match, MsgDcm2nii] = xASL_io_dcm2nii(scanpath, destdir, scan_name, 'DicomFilter', imPar.dcmExtFilter, 'Verbose', imPar.bVerbose, 'Overwrite', imPar.bOverwrite, 'Version', imPar.dcm2nii_version, 'x', x);
								
								% If dcm2nii produced a warning or error, catch this & store it
								if ~isempty(MsgDcm2nii) && ~isempty(regexpi(MsgDcm2nii,'.*(error).*')) % if it contains a warning/error
									dcm2niiCatchedErrors = xASL_bids_CatchErrors('xASL_io_dcm2nii', MsgDcm2nii, dbstack, ['dcm2nii_' imPar.dcm2nii_version], pwd, scan_name, scanpath, destdir, dcm2niiCatchedErrors, imPar);
								end
								
							catch ME
								dcm2niiCatchedErrors = xASL_bids_CatchErrors(ME.identifier, ME.message, [], [], [], scan_name, scanpath, destdir, dcm2niiCatchedErrors, imPar, ME.stack);
								
								if imPar.bVerbose; warning(['dcm2nii ' scanpath ' crashed, skipping']); end
								if imPar.bVerbose; warning('Check whether the scan is complete'); end
								first_match = xASL_adm_GetFileList(scanpath, ['.*' imPar.dcmExtFilter],'FPList',[0 Inf]);
								if  ~isempty(first_match); first_match = first_match{1}; end
							end
						end
						
						%% In case of a single NII ASL file loaded from PAR/REC, we need to shuffle the dynamics from CCCC...LLLL order to CLCLCLCL... order
						[~,~,scanExtension] = xASL_fileparts(scanpath);
						if ~isempty(regexpi(scanExtension, '^\.(par|rec)$')) && length(nii_files)==1 && ~isempty(regexpi(scan_name, 'ASL'))
							% For a PAR/REC files that produces a single ASL4D NIFTI
							imASL = xASL_io_Nifti2Im(nii_files{1});
							% If multiple dynamics
							if size(imASL,4) > 1
								% Then reshuffle them
								imASLreordered = zeros(size(imASL));
								imASLreordered(:,:,:,1:2:end) = imASL(:,:,:,1:ceil(size(imASL,4)/2));
								imASLreordered(:,:,:,2:2:end) = imASL(:,:,:,ceil(size(imASL,4)/2)+1:end);
								xASL_io_SaveNifti(nii_files{1},nii_files{1},imASLreordered);
							end
						end
						% Merge NIfTIs if there are multiples
						% For ASL or M0, merge multiple files
						if length(nii_files)>1
							if ~isempty(strfind(scan_name,'ASL4D'))
								nii_files = xASL_bids_MergeNifti(nii_files, 'ASL');
							elseif  ~isempty(strfind(scan_name,'M0'))
								nii_files = xASL_bids_MergeNifti(nii_files, 'M0');
							end
						end
						
						% Extract relevant parameters from nifti header and append to summary file
						summary_line = xASL_bids_AppendNiftiParameters(nii_files);
						converted_scans(iSubject, iSession, iScan) = 1;
					end
					
					% extract relevant parameters from dicom header, if not
					% already exists
					% Find JSONpath that is there already
					SavePathJSON = {};
					SavePathJSON{1} = fullfile(destdir, [scan_name '.json']);
					SavePathJSON{2} = fullfile(destdir, [session_name '.json']);
					for iPath=1:length(nii_files)
						% now we add the path only if it didnt exist already in this list
						tmpNewPath = [nii_files{iPath}(1:end-4) '.json'];
						if ~max(cellfun(@(y) strcmp(y, tmpNewPath), SavePathJSON))
							SavePathJSON{end+1} = tmpNewPath;
						end
					end
					
					for iPath=1:length(SavePathJSON)
						if exist(SavePathJSON{iPath}, 'file') && ~isempty(first_match)
							[~, ~, fext] = fileparts(first_match);
							if  strcmpi(fext,'.PAR')
								parms = xASL_bids_Par2JSON(first_match, SavePathJSON{iPath});
							elseif strcmpi(fext,'.nii')
								parms = [];
							elseif imPar.bMatchDirectories
								Fpath  = fileparts(first_match);
								[parms, pathDcmDict] = xASL_bids_Dicom2JSON(imPar, Fpath, SavePathJSON{iPath}, imPar.dcmExtFilter, bUseDCMTK, pathDcmDict);
								clear Fpath Ffile Fext
							else
								[parms, pathDcmDict] = xASL_bids_Dicom2JSON(imPar, first_match, SavePathJSON{iPath}, imPar.dcmExtFilter, bUseDCMTK, pathDcmDict);
							end
						end
					end
					
					% correct nifti rescale slope if parms.RescaleSlopeOriginal =~1
					% but nii.dat.scl_slope==1 (this can happen in case of
					% hidden scale slopes in private Philips header,
					% that is dealt with by xASL_bids_Dicom2JSON but not by
					% dcm2niiX
					
					if ~isempty(nii_files) && exist('parms','var')
						[TempLine, PrintDICOMFields] = xASL_bids_AppendParmsParameters(parms);
						summary_line = [summary_line TempLine];
					end
					
					if bClone2Source % make a copy of analysisdir in sourcedir
						if ~isempty(nii_files)
							for iFile=1:length(nii_files)
								% replace 'analysis' by 'source'
								[iStart, iEnd] = regexp(nii_files{iFile}, 'analysis');
								DestPath = [nii_files{iFile}(1:iStart-1) 'source' nii_files{iFile}(iEnd+1:end)];
								xASL_Copy(nii_files{iFile}, DestPath, true);
								% do the same for other extensions
								Extensions = {'.json' '_parms.json'};
								for iExt=1:length(Extensions)
									[Fpath, Ffile] = xASL_fileparts(nii_files{iFile});
									CopyPath = fullfile(Fpath, [Ffile Extensions{iExt}]);
									[Fpath, Ffile] = xASL_fileparts(DestPath);
									DestPath = fullfile(Fpath, [Ffile Extensions{iExt}]);
									if xASL_exist(CopyPath)
										xASL_Copy(CopyPath, DestPath, true);
									end
								end
							end
						end
					end
					
					% Copy single dicom as QC placeholder
					if bCopySingleDicoms && ~isempty(first_match)
						xASL_Copy(first_match, fullfile(destdir, ['DummyDicom_' scan_name '.dcm']), imPar.bOverwrite, imPar.bVerbose);
					end
					
					% store the summary info so it can be sorted and printed below
					summary_lines{iSubject, iVisit, iSession, iScan} = summary_line;
				end % scansIDs
			end % sessionIDs
		end % visitIDs
	end % subjectIDs
	
    % Create summary file
    xASL_bids_CreateSummaryFile(imPar, PrintDICOMFields, ...
                                converted_scans, ...
                                skipped_scans, ...
                                missing_scans, ...
                                subjectIDs, visitIDs, scanNames, ...
                                summary_lines, ...
                                nSubjects, nVisits, nSessions, fid_summary);
    
	% cleanup
	if ~bUseDCMTK || isempty(pathDcmDict)
		dicomdict('factory');
	end
	diary('off');
	
	if ~isempty(fields(dcm2niiCatchedErrors))
		fclose all;
		SavePath = fullfile(imPar.AnalysisRoot, 'dcm2niiCatchedErrors.mat');
		SaveJSON = fullfile(imPar.AnalysisRoot, 'dcm2niiCatchedErrors.json');
		xASL_delete(SavePath);
		xASL_delete(SaveJSON);
		save(SavePath,'dcm2niiCatchedErrors');
		spm_jsonwrite(SaveJSON, dcm2niiCatchedErrors);
	end
	
	fprintf('\n');
end



