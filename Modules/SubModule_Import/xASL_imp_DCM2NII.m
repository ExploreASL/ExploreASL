function xASL_imp_DCM2NII(imPar, x)
%xASL_imp_DCM2NII Run the dcm2nii part of the import.
%
% FORMAT: xASL_imp_DCM2NII(imPar, x)
% 
% INPUT:
%   imPar              - JSON file with structure with import parameters (REQUIRED, STRUCT)
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
% EXAMPLE:     xASL_imp_DCM2NII(imPar, x);
% __________________________________
% Copyright 2015-2021 ExploreASL

    
    %% 1. Initialize defaults of dcm2nii
    dcm2niiCatchedErrors = struct; % initialization
    if x.modules.import.settings.bCheckPermissions
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
	
	if x.modules.import.settings.bCheckPermissions
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
	x.modules.import.pathDcmDict = fullfile(x.MyPath,'External','xASL_DICOMLibrary.txt');
	if ~x.modules.import.settings.bUseDCMTK
		% -----------------------------------------------------------------------------
		% Initialize dicom dictionary by appending private philips stuff to a temporary copy
		% -----------------------------------------------------------------------------
		dicomdict('set', x.modules.import.pathDcmDict);
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
	x.modules.import.listsIDs.vSubjectIDs = tokens(:,imPar.tokenOrdering(1)); % cell vector with extracted subject IDs (for all sessions and scans)
	
	%% VISITS
	if imPar.tokenOrdering(2)==0
		% a zero means: no visits applicable
		x.modules.import.settings.bUseVisits = false;
		x.modules.import.listsIDs.vVisitIDs = cellfun(@(y) '1', x.modules.import.listsIDs.vSubjectIDs, 'UniformOutput', false); % each subject has a single visit
		imPar.tokenVisitAliases = {'^1$', '_1'};
	else
		x.modules.import.settings.bUseVisits = true;
		x.modules.import.listsIDs.vVisitIDs = tokens(:,imPar.tokenOrdering(2)); % cell vector with extracted session IDs (for all subjects and scans)
	end
	
	%% SESSIONS
	if imPar.tokenOrdering(3)==0
		% a zero means: no sessions applicable
		x.modules.import.settings.bUseSessions = false;
		x.modules.import.listsIDs.vSessionIDs = cellfun(@(y) '1', x.modules.import.listsIDs.vSubjectIDs, 'UniformOutput', false); % each subject-visit has a single session
		imPar.tokenSessionAliases = {'^1$', 'ASL_1'};
	else
		x.modules.import.settings.bUseSessions = true;
		x.modules.import.listsIDs.vSessionIDs = tokens(:,imPar.tokenOrdering(3)); % cell vector with extracted session IDs (for all subjects and scans)
	end
	
	%% SCANTYPES
	x.modules.import.listsIDs.vScanIDs = tokens(:,imPar.tokenOrdering(4)); % cell vector with extracted scan IDs (for all subjects and sessions)
	
	% Convert the vectors to unique & sort sets by the output aliases
	x.modules.import.listsIDs.subjectIDs  = sort(unique(x.modules.import.listsIDs.vSubjectIDs));
	x.modules.import.numOf.nSubjects = length(x.modules.import.listsIDs.subjectIDs);
	
	x.modules.import.listsIDs.visitIDs  = unique(x.modules.import.listsIDs.vVisitIDs);
	% sort by output
	if length(x.modules.import.listsIDs.visitIDs)>1
		for iV=1:length(x.modules.import.listsIDs.visitIDs)
			IDrow(iV) = find(cellfun(@(y) strcmp(y,x.modules.import.listsIDs.visitIDs{iV}), imPar.tokenVisitAliases(:,1)));
		end
		x.modules.import.listsIDs.visitIDs = x.modules.import.listsIDs.visitIDs(IDrow);
	end
	x.modules.import.numOf.nVisits = length(x.modules.import.listsIDs.visitIDs);
	
	x.modules.import.listsIDs.sessionIDs  = sort(unique(x.modules.import.listsIDs.vSessionIDs));
	x.modules.import.numOf.nSessions = length(x.modules.import.listsIDs.sessionIDs);
	
	x.modules.import.listsIDs.scanIDs = sort(unique(lower(x.modules.import.listsIDs.vScanIDs)));
	x.modules.import.numOf.nScans = length(x.modules.import.listsIDs.scanIDs);
	
	%% VISIT NAMES
	if isempty(imPar.visitNames)
		if isempty(x.modules.import.listsIDs.visitIDs)
			imPar.visitNames = cell(x.modules.import.numOf.nVisits,1);
			for kk=1:x.modules.import.numOf.nVisits
				imPar.visitNames{kk}=sprintf('ASL_%g',kk);
			end
		else
			imPar.visitNames = x.modules.import.listsIDs.visitIDs;
		end
	end
	
	%% SESSION NAMES
	% optionaly we can have human readble session names; by default they are the same as the original tokens in the path
	if isempty(imPar.sessionNames)
		if isempty(x.modules.import.listsIDs.sessionIDs)
			imPar.sessionNames = cell(x.modules.import.numOf.nSessions,1);
			for kk=1:x.modules.import.numOf.nSessions
				imPar.sessionNames{kk}=sprintf('ASL_%g',kk);
			end
		else
			imPar.sessionNames = x.modules.import.listsIDs.sessionIDs;
		end
	end
	
	%% SCAN NAMES
	x.modules.import.scanNames = x.modules.import.listsIDs.scanIDs;
	
	%% 6. Sanity check for missing elements
	if x.modules.import.numOf.nSubjects==0
		error('No subjects')
	end
	if x.modules.import.numOf.nVisits==0
		error('No visits')
	end
	if x.modules.import.numOf.nSessions==0
		error('No sessions')
	end
	if x.modules.import.numOf.nScans==0
		error('No scans')
	end
	
	% preallocate space for (global) counts
	x.modules.import.globalCounts.converted_scans = zeros(x.modules.import.numOf.nSubjects,...
                                                          x.modules.import.numOf.nVisits,...
                                                          x.modules.import.numOf.nSessions,...
                                                          x.modules.import.numOf.nScans,'uint8'); % keep a count of all individual scans
	x.modules.import.globalCounts.skipped_scans = zeros(x.modules.import.numOf.nSubjects,...
                                                        x.modules.import.numOf.nVisits,...
                                                        x.modules.import.numOf.nSessions,...
                                                        x.modules.import.numOf.nScans,'uint8'); % keep a count of all individual scans
	x.modules.import.globalCounts.missing_scans = zeros(x.modules.import.numOf.nSubjects,...
                                                        x.modules.import.numOf.nVisits,...
                                                        x.modules.import.numOf.nSessions,...
                                                        x.modules.import.numOf.nScans,'uint8'); % keep a count of all individual scans
	
	% define a cell array for storing info for parameter summary file
	xASL_adm_CreateDir(imPar.AnalysisRoot);
	x.modules.import.summary_lines = cell(x.modules.import.numOf.nSubjects,...
                                          x.modules.import.numOf.nVisits,...
                                          x.modules.import.numOf.nSessions,...
                                          x.modules.import.numOf.nScans);
	
    
	%% 7. Import subject by subject, visit by visit, session by session, scan by scan
	fprintf('%s\n', 'Running import (i.e. dcm2niiX)');
    
    % Iterate over subjects
    for iSubject=1:x.modules.import.numOf.nSubjects
        [imPar, x.modules.import.summary_lines, PrintDICOMFields, x.modules.import.globalCounts, x.modules.import.scanNames, dcm2niiCatchedErrors, x.modules.import.pathDcmDict] = ...
            xASL_imp_DCM2NII_Subject(x, imPar, iSubject, matches, dcm2niiCatchedErrors);
    end
	
    % Create summary file
    xASL_imp_CreateSummaryFile(imPar, PrintDICOMFields, x, fid_summary);
    
	% cleanup
	if ~x.modules.import.settings.bUseDCMTK || isempty(x.modules.import.pathDcmDict)
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



