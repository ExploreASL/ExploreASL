function xASL_imp_DCM2NII(imPar, bCopySingleDicoms, bUseDCMTK, bCheckPermissions, bClone2Source, x)
%xASL_imp_DCM2NII Run the dcm2nii part of the import.
%
% FORMAT: xASL_imp_DCM2NII(imPar, bCopySingleDicoms, bUseDCMTK, bCheckPermissions, bClone2Source,x)
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
%
% 1. Initialize defaults of dcm2nii
% 2. Create the basic folder structure for sourcedata & derivative data
% 3. Here we try to fix backwards compatibility, but this may break
% 4. Redirect output to a log file
% 5. Start with defining the subjects, visits, sessions (i.e. BIDS runs) and scans (i.e. ScanTypes) by listing or typing
% 6. Sanity check for missing elements
% 7. Import subject by subject, visit by visit, session by session, scan by scan
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     xASL_imp_DCM2NII(imPar, bCopySingleDicoms, bUseDCMTK, bCheckPermissions, bClone2Source,x);
% __________________________________
% Copyright 2015-2021 ExploreASL

    
    %% 1. Initialize defaults of dcm2nii
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
	
	%% 2. Create the basic folder structure for sourcedata & derivative data
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
	
	%% 3. Here we try to fix backwards compatibility, but this may break
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
	
	%% 4. Redirect output to a log file
	diary_filepath = fullfile(imPar.AnalysisRoot, ['import_log_' imPar.studyID '_' datestr(now,'yyyymmdd_HHMMSS') '.txt']);
	diary(diary_filepath);
	
	if imPar.bMatchDirectories
		strLookFor = 'Directories';
	else
		strLookFor = 'Files';
	end
	
	
	%% 5. Start with defining the subjects, visits, sessions (i.e. BIDS runs) and scans (i.e. ScanTypes) by listing or typing
	
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
	
	%% VISITS
	if imPar.tokenOrdering(2)==0
		% a zero means: no visits applicable
		bUseVisits = false;
		vVisitIDs = cellfun(@(y) '1', vSubjectIDs, 'UniformOutput', false); % each subject has a single visit
		imPar.tokenVisitAliases = {'^1$', '_1'};
	else
		bUseVisits = true;
		vVisitIDs = tokens(:,imPar.tokenOrdering(2)); % cell vector with extracted session IDs (for all subjects and scans)
	end
	
	%% SESSIONS
	if imPar.tokenOrdering(3)==0
		% a zero means: no sessions applicable
		bUseSessions = false;
		vSessionIDs = cellfun(@(y) '1', vSubjectIDs, 'UniformOutput', false); % each subject-visit has a single session
		imPar.tokenSessionAliases = {'^1$', 'ASL_1'};
	else
		bUseSessions = true;
		vSessionIDs = tokens(:,imPar.tokenOrdering(3)); % cell vector with extracted session IDs (for all subjects and scans)
	end
	
	%% SCANTYPES
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
	
	%% VISIT NAMES
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
	
	%% SESSION NAMES
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
	
	%% SCAN NAMES
	scanNames = scanIDs;
	
	%% 6. Sanity check for missing elements
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
	globalCounts.converted_scans = zeros(nSubjects,nVisits,nSessions,nScans,'uint8'); % keep a count of all individual scans
	globalCounts.skipped_scans = zeros(nSubjects,nVisits,nSessions,nScans,'uint8'); % keep a count of all individual scans
	globalCounts.missing_scans = zeros(nSubjects,nVisits,nSessions,nScans,'uint8'); % keep a count of all individual scans
	
	% define a cell array for storing info for parameter summary file
	xASL_adm_CreateDir(imPar.AnalysisRoot);
	summary_lines = cell(nSubjects,nVisits,nSessions,nScans);
	
    
	%% 7. Import subject by subject, visit by visit, session by session, scan by scan
	fprintf('%s\n', 'Running import (i.e. dcm2niiX)');
    
    % Create ID struct
    listsIDs.vSubjectIDs = vSubjectIDs;
    listsIDs.vVisitIDs = vVisitIDs;
    listsIDs.vSessionIDs = vSessionIDs;
    listsIDs.vScanIDs = vScanIDs;
    listsIDs.subjectIDs = subjectIDs;
    listsIDs.visitIDs = visitIDs;
    listsIDs.sessionIDs = sessionIDs;
    listsIDs.scanIDs = scanIDs;
    
    % Create number of struct
    numOf.nSubjects = nSubjects;
    numOf.nVisits = nVisits;
    numOf.nSessions = nSessions;
    numOf.nScans = nScans;
    
    % Create settings struct
    settings.bUseVisits = bUseVisits;
    settings.bClone2Source = bClone2Source;
    settings.bUseDCMTK = bUseDCMTK;
    settings.bCopySingleDicoms = bCopySingleDicoms;
    
    % Iterate over subjects
    for iSubject=1:nSubjects
        [imPar, summary_lines, PrintDICOMFields, globalCounts, dcm2niiCatchedErrors, pathDcmDict] = ...
            xASL_imp_DCM2NII_Subject(x, imPar, listsIDs, numOf, settings, globalCounts, iSubject, summary_lines, matches, dcm2niiCatchedErrors, pathDcmDict);
    end
	
    % Create summary file
    xASL_imp_CreateSummaryFile(imPar, numOf, listsIDs, PrintDICOMFields, globalCounts, scanNames, summary_lines, fid_summary);
    
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



