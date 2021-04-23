function [imPar, summary_lines, PrintDICOMFields, globalCounts, dcm2niiCatchedErrors, pathDcmDict] = xASL_imp_DCM2NII_Subject(x, imPar, listsIDs, numOf, settings, globalCounts, iSubject, summary_lines, matches, dcm2niiCatchedErrors, pathDcmDict)
%xASL_imp_DCM2NII_Subject Run DCM2NII for one individual subject.
%
% FORMAT: [imPar, summary_lines, PrintDICOMFields, globalCounts, dcm2niiCatchedErrors, pathDcmDict] = xASL_imp_DCM2NII_Subject(x, imPar, listsIDs, numOf, settings, globalCounts, iSubject, summary_lines, matches, dcm2niiCatchedErrors, pathDcmDict)
% 
% INPUT:
%   x                      - ExploreASL x structure (REQUIRED, STRUCT)
%   imPar                  - Structure with import parameters (REQUIRED, STRUCT)
%   listsIDs               - Lists of IDs (REQUIRED, STRUCT)
%   numOf                  - Number of visits, sessions, scans etc. (REQUIRED, STRUCT)
%   settings               - Boolean settings (REQUIRED, STRUCT)
%   globalCounts           - Converted, skipped & missing scans (REQUIRED, STRUCT)
%   iSubject               - Current subject (REQUIRED, INTEGER)
%   summary_lines          - Summary lines (REQUIRED, CELL ARRAY)
%   matches                - Matches (REQUIRED, CELL ARRAY)
%   dcm2niiCatchedErrors   - DCM2NII catched errors (REQUIRED, STRUCT)
%   pathDcmDict            - Path to DCM dictionary (REQUIRED, CHAR ARRAY)
%
% OUTPUT:
%   imPar                  - Structure with import parameters 
%   summary_lines          - Summary lines
%   PrintDICOMFields       - Print DICOM fields
%   globalCounts           - Converted, skipped & missing scans
%   dcm2niiCatchedErrors   - DCM2NII catched errors
%   pathDcmDict            - Path to DCM dictionary
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Run DCM2NII for one individual subject.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     [imPar, summary_lines, PrintDICOMFields, globalCounts, dcm2niiCatchedErrors, pathDcmDict] = ...
%               xASL_imp_DCM2NII_Subject(x, imPar, listsIDs, numOf, settings, globalCounts, iSubject, summary_lines, matches, dcm2niiCatchedErrors, pathDcmDict)
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% Run DCM2NII for one individual subject
    
    separatorline = repmat(char('+'),1,80);
    subjectID = listsIDs.subjectIDs{iSubject};

    for iVisit=1:numOf.nVisits
        visitID = listsIDs.visitIDs{iVisit};

        % convert visit ID to a suitable name
        if size(imPar.tokenVisitAliases,2)==2
            iAlias = find(~cellfun(@isempty,regexp(listsIDs.visitIDs{iVisit},imPar.tokenVisitAliases(:,1),'once')));
            if ~isempty(iAlias)
                imPar.visitNames{iVisit} = imPar.tokenVisitAliases{iAlias,2};
            end
        end

        if settings.bUseVisits % only pad VisitID _1 _2 _3 etc if there are visits specified
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
        for iSession=1:numOf.nSessions
            sessionID = listsIDs.sessionIDs{iSession};

            % convert session ID to a suitable name
            if size(imPar.tokenSessionAliases,2)==2
                iAlias = find(~cellfun(@isempty,regexp(sessionID,imPar.tokenSessionAliases(:,1),'once')));
                if ~isempty(iAlias)
                    imPar.sessionNames{iSession} = imPar.tokenSessionAliases{iAlias,2};
                end
            end

            for iScan=1:numOf.nScans
                scanID = listsIDs.scanIDs{iScan};
                summary_line = [];
                first_match = [];
                summary_lines{iSubject,iVisit,iSession,iScan} = 'n/a';

                if ~imPar.bVerbose % if not verbose, track % progress
                    CounterN = (iSession-1)*numOf.nScans+iScan;
                    CounterT = numOf.nSessions*numOf.nScans;
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
                        dcm2niiCatchedErrors = xASL_imp_CatchErrors('isempty(iAlias)', WarningMessage, dbstack, mfilename, pwd, scan_name, scanpath, destdir, dcm2niiCatchedErrors, imPar);
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
                iMatch = find(strcmp(listsIDs.vSubjectIDs,subjectID) & strcmp(listsIDs.vVisitIDs, xASL_adm_CorrectName(visitID,2,'_')) & strcmp(listsIDs.vSessionIDs,sessionID) & strcmpi(listsIDs.vScanIDs,scanID) ); % only get the matching session
                if isempty(iMatch)
                    % only report as missing if we need a scan for each session (i.e. ASL)
                    if sum(globalCounts.converted_scans(iSubject,iVisit,:,iScan))==0
                        WarningMessage = ['Missing scan: ' [subjectID imPar.visitNames{iVisit}] ', ' num2str(iSession) ', ' scan_name];
                        if imPar.bVerbose; warning(WarningMessage); end
                        globalCounts.missing_scans(iSubject, iVisit, iSession, iScan) = 1;
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

                if bOneScanIsEnough && sum(globalCounts.converted_scans(iSubject,iVisit,:,iScan))~=0
                    % one scan is enough, so skip this one if there was already a scan converted of this type (i.e. T1)
                    if imPar.bVerbose; fprintf('Skipping scan: %s, %s, %s\n',[subjectID imPar.visitNames{iVisit}],session_name,scan_name); end
                    bSkipThisOne = true;
                    destdir = []; % just in case
                end

                % start the conversion if this scan should not be skipped
                if bSkipThisOne
                    summary_line = sprintf(',"skipped",,,,,,,,');
                    globalCounts.skipped_scans(iSubject, iVisit, iSession, iScan) = 1;
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
                                dcm2niiCatchedErrors = xASL_imp_CatchErrors('xASL_io_dcm2nii', MsgDcm2nii, dbstack, ['dcm2nii_' imPar.dcm2nii_version], pwd, scan_name, scanpath, destdir, dcm2niiCatchedErrors, imPar);
                            end

                        catch ME
                            dcm2niiCatchedErrors = xASL_imp_CatchErrors(ME.identifier, ME.message, [], [], [], scan_name, scanpath, destdir, dcm2niiCatchedErrors, imPar, ME.stack);

                            if imPar.bVerbose; warning(['dcm2nii ' scanpath ' crashed, skipping']); end
                            if imPar.bVerbose; warning('Check whether the scan is complete'); end
                            first_match = xASL_adm_GetFileList(scanpath, ['.*' imPar.dcmExtFilter],'FPList',[0 Inf]);
                            if  ~isempty(first_match); first_match = first_match{1}; end
                        end
                    end

                    %% In case of a single NII ASL file loaded from PAR/REC, we need to shuffle the dynamics from CCCC...LLLL order to CLCLCLCL... order
                    [nii_files, summary_line, globalCounts, ASLContext] = xASL_imp_DCM2NII_Subject_ShuffleTheDynamics(globalCounts, scanpath, scan_name, nii_files, iSubject, iSession, iScan);
                    
                    
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
                
                % Store JSON files
                [parms, pathDcmDict] = xASL_imp_DCM2NII_Subject_StoreJSON(imPar, SavePathJSON, first_match, settings.bUseDCMTK, pathDcmDict);
                
                % Get real ASLContext
                realASLContext = xASL_imp_DCM2NII_Subject_RealASLContext(parms);
                
                % Check Context
                if ~strcmp(ASLContext,realASLContext)
                    fprintf('Context mismatch, will switch NIFTI ordering...\n');
                    xASL_imp_DCM2NII_Subject_SwitchNIIs(nii_files);
                end

                % correct nifti rescale slope if parms.RescaleSlopeOriginal =~1
                % but nii.dat.scl_slope==1 (this can happen in case of
                % hidden scale slopes in private Philips header,
                % that is dealt with by xASL_bids_Dicom2JSON but not by
                % dcm2niiX

                if ~isempty(nii_files) && exist('parms','var')
                    [TempLine, PrintDICOMFields] = xASL_imp_AppendParmsParameters(parms);
                    summary_line = [summary_line TempLine];
                else
                    PrintDICOMFields = [];
                end

                % Make a copy of analysisdir in sourcedir
                xASL_imp_DCM2NII_Subject_CopyAnalysisDir(nii_files, settings.bClone2Source)

                % Copy single dicom as QC placeholder
                if settings.bCopySingleDicoms && ~isempty(first_match)
                    xASL_Copy(first_match, fullfile(destdir, ['DummyDicom_' scan_name '.dcm']), imPar.bOverwrite, imPar.bVerbose);
                end

                % store the summary info so it can be sorted and printed below
                summary_lines{iSubject, iVisit, iSession, iScan} = summary_line;
            end % scansIDs
        end % sessionIDs
    end % visitIDs

        
end


%% Store JSON files
function [parms, pathDcmDict] = xASL_imp_DCM2NII_Subject_StoreJSON(imPar, SavePathJSON, first_match, bUseDCMTK, pathDcmDict)

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

end

%% Make a copy of analysisdir in sourcedir
function xASL_imp_DCM2NII_Subject_CopyAnalysisDir(nii_files, bClone2Source)

    if bClone2Source
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

end


%% Shuffle the dynamics
function [nii_files, summary_line, globalCounts, ASLContext] = xASL_imp_DCM2NII_Subject_ShuffleTheDynamics(globalCounts, scanpath, scan_name, nii_files, iSubject, iSession, iScan)

    % Set to true if you want to print information for debug purposes
    bVerbose = true;
    
    % Print NIFTI files
    if bVerbose
        fprintf('Shuffle the dynamics...\n');
        for iNifti = 1:size(nii_files,2)
            fprintf('%s\n',nii_files{1,iNifti});
        end
    end

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
            xASL_io_SaveNifti(nii_files{1}, nii_files{1}, imASLreordered);
        end
    end
    
    % Merge NIfTIs if there are multiples
    % For ASL or M0, merge multiple files
    if length(nii_files)>1
        if ~isempty(strfind(scan_name,'ASL4D'))
            [nii_files,ASLContext] = xASL_bids_MergeNifti(nii_files, 'ASL');
        elseif  ~isempty(strfind(scan_name,'M0'))
            [nii_files,ASLContext] = xASL_bids_MergeNifti(nii_files, 'M0');
        end
    end    
    
    % Extract relevant parameters from nifti header and append to summary file
    summary_line = xASL_imp_AppendNiftiParameters(nii_files);
    globalCounts.converted_scans(iSubject, iSession, iScan) = 1;

end


%% Get real ASL Context
function realASLContext = xASL_imp_DCM2NII_Subject_RealASLContext(parms)

    % Fallback
    realASLContext = '';

    if isfield(parms,'InstanceNumber') && isfield(parms,'ImageTypes')
        if size(parms.InstanceNumber,2)==size(parms.ImageTypes,2)
            % Get table
            for iNum = 1:size(parms.InstanceNumber,2)
                tableInNumImType{iNum,1} = parms.InstanceNumber(1,iNum);
                tableInNumImType{iNum,2} = parms.ImageTypes(1,iNum);
            end
            % Sort rows
            tableInNumImType = sortrows(tableInNumImType,1);
            % Get ASL context
            realASLContextCell = unique([tableInNumImType{:,2}]);
            for iContext = 1:size(realASLContextCell,2)
                parStruct.ImageType = realASLContextCell{1,iContext};
                thisContext = xASL_bids_determineImageTypeGE(parStruct);
                realASLContext = [realASLContext ',' thisContext];
            end
            realASLContext = realASLContext(2:end);
        end
    end


end


%% Get real ASL Context
function xASL_imp_DCM2NII_Subject_SwitchNIIs(nii_files)

    % Load image
    imASL = xASL_io_Nifti2Im(nii_files{1});

    % Then reshuffle them
    imASLreordered = zeros(size(imASL));
    imASLreordered(:,:,:,1:2:end) = imASL(:,:,:,1:ceil(size(imASL,4)/2));
    imASLreordered(:,:,:,2:2:end) = imASL(:,:,:,ceil(size(imASL,4)/2)+1:end);

    % Save image
    xASL_io_SaveNifti(nii_files{1}, nii_files{1}, imASLreordered);

end


