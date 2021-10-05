function xASL_imp_DCM2NII(imPar, x, newLogging)
%xASL_imp_DCM2NII Run the dcm2nii part of the import.
%
% FORMAT: xASL_imp_DCM2NII(imPar, x)
% 
% INPUT:
%   imPar              - JSON file with structure with import parameters (REQUIRED, STRUCT)
%   x                  - ExploreASL x structure (REQUIRED, STRUCT)
%   newLogging         - Option to create new xASL_module_Import log (BOOLEAN, OPTION, DEFAULT = true)
%                        (using the false option you will get the well known dcm2niix log)
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

    %% 0. Input check
    if nargin<3 || isempty(newLogging)
        newLogging = true;
    end
    
    %% 1. Initialize defaults of dcm2nii
    fprintf('================================== DICOM to NIFTI CONVERSION =================================\n');
    
    % Create the temp directory for DCM2NII
    xASL_adm_CreateDir(imPar.TempRoot);
    
    % Initialization of an empty catched errors struct
    dcm2niiCatchedErrors = struct;
    
    % Basic import checks before execution
    [x,imPar] = xASL_imp_CheckImportSettings(x,imPar);
	
	%% 2. Create the basic folder structure for sourcedata & derivative data
    if ~exist(imPar.RawRoot, 'dir')
        warning(['Could not find ' imPar.RawRoot ', trying to find a different folder instead...']);
        
        % find any folder except for temp, sourcedata, rawdata, derivatives
        % xASL_adm_GetFileList uses regular expressions, to create a nice list of foldernames,
        % with/without FullPath (FPList), with/without recursive (FPListRec)
        % very powerful once you know how these work
        FolderNames = xASL_adm_GetFileList(fullfile(imPar.RawRoot, imPar.studyID), ...
            '^(?!(temp|derivatives|source|sourcedata)).*$', 'FPList', [0 Inf], true);
        
        if length(FolderNames)==1
            imPar.RawRoot = FolderNames{1};
            fprintf('%s\n', ['Found ' imPar.RawRoot ' as sourcedata folder']);
        else
            error('Couldnt find a sourcedata folder, please rename one, or move other folders');
        end
    end
	
    % Check access rights of temp and rawdata directories
	if x.modules.import.settings.bCheckPermissions
		xASL_adm_CheckPermissions(imPar.RawRoot, false); % don't need execution permisions
		xASL_adm_CheckPermissions(imPar.TempRoot, false);  % don't need execution permisions
	end
	
	%% 3. Here we try to fix backwards compatibility
	if length(imPar.tokenOrdering)==3
        % Backwards compatibility visits
		if size(imPar.tokenOrdering,1) > 1
			% Vertical vector
			imPar.tokenOrdering = imPar.tokenOrdering';
		end
		imPar.tokenOrdering = [imPar.tokenOrdering(1) 0 imPar.tokenOrdering(2:3)]; % insert Visits(none)
	elseif length(imPar.tokenOrdering)==2 
        % Backwards compatibility Visits & sessions
		imPar.tokenOrdering = [imPar.tokenOrdering(1) 0 0 imPar.tokenOrdering(2)]; % insert sessions & visits(none)
	end
	
	% Path to the dictionary to initialize - we need to keep track if the dictionary has been set, because
	% Dicominfo can be used despite bUSEDCMTK==1 when DCMTK fails
    x.modules.import.pathDcmDict = fullfile(x.opts.MyPath,'External','xASL_DICOMLibrary.txt');
    if ~x.modules.import.settings.bUseDCMTK
        % Initialize dicom dictionary by appending private philips stuff to a temporary copy
        dicomdict('set', x.modules.import.pathDcmDict);
    end
	
    % Initialize to be able to catch errors and close if valid
	fid_summary = -1;
	
	% change dcmnii_version for PARREC if needed
	if ~isempty(strfind(char(imPar.folderHierarchy(end)),'PAR'))
		imPar.dcm2nii_version = '20101105';
	end
	
	%% 4. Redirect output to a log file
    if newLogging
        diary_filepath = fullfile(x.opts.DatasetRoot,'xASL_module_Import.log');
    else
        diary_filepath = fullfile(imPar.TempRoot, ['import_log_' imPar.studyID '_' datestr(now,'yyyymmdd_HHMMSS') '.txt']);
    end
	diary(diary_filepath);
	
    % Check if directories are files are supposed to be matched
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
    
    % Print matching files
	if isempty(matches)
		warning('No matching files, skipping');
		return;
	elseif imPar.bVerbose
		fprintf('\nMatching files (#=%g):\n',length(matches));
        for iMatch=1:size(matches,1)
            fprintf('%s\n', matches{iMatch,1});
        end
	end
	
	% Copy the columns into named vectors. This construction allows for arbitrary directory hierarchies.
	% Make sure to select the columns that correspond to the folder ordering defined using the regular expressions above.
	
    % Define Subjects (so subjects may be repeated here)
    % vSubjectIDs: cell vector with extracted subject IDs (for all visits, sessions and scans; 
	x.modules.import.listsIDs.vSubjectIDs = tokens(:,imPar.tokenOrdering(1));
	
    
    %% Determine structure from sourcedata
    [x,imPar] = xASL_imp_DetermineStructureFromSourcedata(x,imPar,tokens);
    
	
	% 6. Sanity check for missing elements
    xASL_imp_DCM2NII_SanityChecks(x);
	
    
    % Preallocate space for (global) counts
    x = xASL_imp_PreallocateGlobalCounts(x);
	
    
	%% 7. Import subject by subject, visit by visit, session by session, scan by scan
	fprintf('\nRunning DCM2NIIX...\n');
    
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
		SavePath = fullfile(imPar.TempRoot, 'dcm2niiCatchedErrors.mat');
		SaveJSON = fullfile(imPar.TempRoot, 'dcm2niiCatchedErrors.json');
		xASL_delete(SavePath);
		xASL_delete(SaveJSON);
		save(SavePath,'dcm2niiCatchedErrors');
		spm_jsonwrite(SaveJSON, dcm2niiCatchedErrors);
	end
	
	fprintf('\n');
    
    
end


%% Basic import checks before execution
function [x,imPar] = xASL_imp_CheckImportSettings(x,imPar)

    if x.modules.import.settings.bCheckPermissions
        dcm2niiDir = fullfile(x.opts.MyPath, 'External', 'MRIcron');
        xASL_adm_CheckPermissions(dcm2niiDir, true); % dcm2nii needs to be executable
    end
    if ~isfield(imPar,'dcm2nii_version') || isempty(imPar.dcm2nii_version)
        % OR for PARREC imPar.dcm2nii_version = '20101105'; THIS IS AUTOMATED BELOW
        imPar.dcm2nii_version = '20190902';
    end
    if ~isfield(imPar,'dcmExtFilter') || isempty(imPar.dcmExtFilter)
        % dcmExtFilter: the last one is because some convertors save files without extension, 
        % but there would be a dot/period before a bunch of numbers
        imPar.dcmExtFilter = '^(.*\.dcm|.*\.img|.*\.IMA|[^.]+|.*\.\d*)$';
    end
	
	if ~isfield(imPar,'SkipSubjectIfExists') || isempty(imPar.SkipSubjectIfExists)
		% allows to skip existing subject folders in the temp folder, when this is set to true,
		% avoiding partly re-importing/converting dcm2niiX when processing has been partly done
		imPar.SkipSubjectIfExists = false;
	else
		warning('Skipping existing subjects in temp folder...');
		fprintf('If you want to overwrite, first remove the full subject folder...');
	end

end


%% Determine structure from sourcedata
function [x,imPar] = xASL_imp_DetermineStructureFromSourcedata(x,imPar,tokens)


    %% VISITS
	if imPar.tokenOrdering(2)==0
		% a zero means: no visits applicable
		x.modules.import.settings.bUseVisits = false;
        % vVisitIDs: each subject has a single visit
		x.modules.import.listsIDs.vVisitIDs = cellfun(@(y) '1', x.modules.import.listsIDs.vSubjectIDs, 'UniformOutput', false);
		imPar.tokenVisitAliases = {'^1$', '_1'};
	else
		x.modules.import.settings.bUseVisits = true;
        % vVisitIDs: cell vector with extracted session IDs (for all subjects, sessions and scans)
		x.modules.import.listsIDs.vVisitIDs = tokens(:,imPar.tokenOrdering(2));
	end
	
	%% SESSIONS
	if imPar.tokenOrdering(3)==0
		% a zero means: no sessions applicable
		x.modules.import.settings.bUseSessions = false;
        % vSessionIDs: each subject-visit has a single session
		x.modules.import.listsIDs.vSessionIDs = cellfun(@(y) '1', x.modules.import.listsIDs.vSubjectIDs, 'UniformOutput', false);
		imPar.tokenSessionAliases = {'^1$', 'ASL_1'};
	else
		x.modules.import.settings.bUseSessions = true;
        % vSessionIDs: Cell vector with extracted session IDs (for all subjects and scans)
		x.modules.import.listsIDs.vSessionIDs = tokens(:,imPar.tokenOrdering(3));
	end
	
	%% SCANTYPES
    
    % vScanIDs: cell vector with extracted scan IDs (for all subjects, visits and sessions)
	x.modules.import.listsIDs.vScanIDs = tokens(:,imPar.tokenOrdering(4)); 
	
	% Convert the vectors to unique & sort sets by the output aliases
	x.modules.import.listsIDs.subjectIDs  = sort(unique(x.modules.import.listsIDs.vSubjectIDs));
	x.modules.import.numOf.nSubjects = length(x.modules.import.listsIDs.subjectIDs);
	x.modules.import.listsIDs.visitIDs  = unique(x.modules.import.listsIDs.vVisitIDs);
    
	% Sort by output
    if length(x.modules.import.listsIDs.visitIDs)>1
        for iV=1:length(x.modules.import.listsIDs.visitIDs)
            IDrow(iV) = find(cellfun(@(y) strcmp(y,x.modules.import.listsIDs.visitIDs{iV}), imPar.tokenVisitAliases(:,1)));
        end
        x.modules.import.listsIDs.visitIDs = x.modules.import.listsIDs.visitIDs(IDrow);
    end
    
    % Get number of visits, session IDs, number of sessions, scan IDs, number of scans
	x.modules.import.numOf.nVisits = length(x.modules.import.listsIDs.visitIDs);
	x.modules.import.listsIDs.sessionIDs  = sort(unique(x.modules.import.listsIDs.vSessionIDs));
	x.modules.import.numOf.nSessions = length(x.modules.import.listsIDs.sessionIDs);
	x.modules.import.listsIDs.scanIDs = sort(unique(lower(x.modules.import.listsIDs.vScanIDs)));
	x.modules.import.numOf.nScans = length(x.modules.import.listsIDs.scanIDs);
	
	%% VISIT NAMES (== session in BIDS)
	if isempty(imPar.visitNames)
		if isempty(x.modules.import.listsIDs.visitIDs)
			imPar.visitNames = cell(x.modules.import.numOf.nVisits,1);
			for kk=1:x.modules.import.numOf.nVisits
				imPar.visitNames{kk} = sprintf('ASL_%g', kk);
			end
		else
			imPar.visitNames = x.modules.import.listsIDs.visitIDs;
		end
	end
	
	%% SESSION NAMES (== run in BIDS)
	% optionally we can have human readble session names; by default they are the same as the original tokens in the path
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


end



%% Sanity check for missing elements
function xASL_imp_DCM2NII_SanityChecks(x)

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

end



%% Preallocate space for (global) counts
function x = xASL_imp_PreallocateGlobalCounts(x)

    % keep a count of all individual scans
    x.modules.import.globalCounts.converted_scans = ...
        zeros(x.modules.import.numOf.nSubjects, x.modules.import.numOf.nVisits, x.modules.import.numOf.nSessions, x.modules.import.numOf.nScans,'uint8');
    % keep a count of all individual scans
    x.modules.import.globalCounts.skipped_scans = ...
        zeros(x.modules.import.numOf.nSubjects, x.modules.import.numOf.nVisits, x.modules.import.numOf.nSessions, x.modules.import.numOf.nScans,'uint8');
    % keep a count of all individual scans
    x.modules.import.globalCounts.missing_scans = ...
        zeros(x.modules.import.numOf.nSubjects, x.modules.import.numOf.nVisits, x.modules.import.numOf.nSessions, x.modules.import.numOf.nScans,'uint8');
	
	% define a cell array for storing info for parameter summary file
	x.modules.import.summary_lines = ...
        cell(x.modules.import.numOf.nSubjects, x.modules.import.numOf.nVisits, x.modules.import.numOf.nSessions, x.modules.import.numOf.nScans);


end



