function [imPar, summary_lines, PrintDICOMFields, globalCounts, scanNames, dcm2niiCatchedErrors, pathDcmDict] = xASL_imp_DCM2NII_Subject(x, imPar, iSubject, matches, dcm2niiCatchedErrors)
%xASL_imp_DCM2NII_Subject Run DCM2NII for one individual subject.
%
% FORMAT: [imPar, summary_lines, PrintDICOMFields, globalCounts, scanNames, dcm2niiCatchedErrors, pathDcmDict] = xASL_imp_DCM2NII_Subject(x, imPar, iSubject, matches, dcm2niiCatchedErrors)
% 
% INPUT:
%   x                      - ExploreASL x structure (REQUIRED, STRUCT)
%   imPar                  - Structure with import parameters (REQUIRED, STRUCT)
%   iSubject               - Current subject (REQUIRED, INTEGER)
%   matches                - Matches (REQUIRED, CELL ARRAY)
%   dcm2niiCatchedErrors   - DCM2NII catched errors (REQUIRED, STRUCT)
%
% OUTPUT:
%   imPar                  - Structure with import parameters 
%   summary_lines          - Summary lines
%   PrintDICOMFields       - Print DICOM fields
%   globalCounts           - Converted, skipped & missing scans
%   scanNames              - Scan names
%   dcm2niiCatchedErrors   - DCM2NII catched errors
%   pathDcmDict            - Path to DCM dictionary
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Run DCM2NII for one individual subject.
%
% 1. Run DCM2NII for one individual subject
% 2. Iterate over visits
% 3. Loop through all sessions
% 4. Iterate over scans
% - 1. Initialize variables (scanID, summary_line, first_match)
% - 2. Convert scan ID to a suitable name and set scan-specific parameters
% - 3. Minimalistic feedback of where we are
% - 4. Now pick the matching one from the folder list
% - 5. Determine input and output paths
% - 6. Start the conversion if this scan should not be skipped
% - 7. Store JSON files
% - 8. In case of a single NII ASL file loaded from PAR/REC, we need to shuffle the dynamics from CCCC...LLLL order to CLCLCLCL... order
% - 9. Make a copy of analysisdir in sourcedir
% - 10. Store the summary info so it can be sorted and printed below
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     [imPar, summary_lines, PrintDICOMFields, globalCounts, scanNames, dcm2niiCatchedErrors, pathDcmDict] = xASL_imp_DCM2NII_Subject(x, imPar, iSubject, matches, dcm2niiCatchedErrors);
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% 1. Run DCM2NII for one individual subject
    
    separatorline = '==============================================================================================';
    summary_lines = x.modules.import.summary_lines;
    pathDcmDict = x.modules.import.pathDcmDict;
    scanNames = x.modules.import.scanNames;
    globalCounts = x.modules.import.globalCounts;
    numOf = x.modules.import.numOf;
    settings = x.modules.import.settings;
    listsIDs = x.modules.import.listsIDs;
    subjectID = listsIDs.subjectIDs{iSubject};
    
    %% 2. Iterate over visits
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

        %% 3. Loop through all sessions
        for iSession=1:numOf.nSessions
            sessionID = listsIDs.sessionIDs{iSession};

            % convert session ID to a suitable name
            if size(imPar.tokenSessionAliases,2)==2
                iAlias = find(~cellfun(@isempty,regexp(sessionID,imPar.tokenSessionAliases(:,1),'once')));
                if ~isempty(iAlias)
                    imPar.sessionNames{iSession} = imPar.tokenSessionAliases{iAlias,2};
                end
            end
            
            %% 4. Iterate over scans
            for iScan=1:numOf.nScans
                
                %% 4.1 Initialize variables (scanID, summary_line, first_match)
                scanID = listsIDs.scanIDs{iScan};
                summary_line = [];
                first_match = [];
                summary_lines{iSubject,iVisit,iSession,iScan} = 'n/a';

                if ~imPar.bVerbose % if not verbose, track % progress
                    CounterN = (iSession-1)*numOf.nScans+iScan;
                    CounterT = numOf.nSessions*numOf.nScans;
                    xASL_TrackProgress(CounterN, CounterT);
                end

                %% 4.2 Convert scan ID to a suitable name and set scan-specific parameters
                if size(imPar.tokenScanAliases,2)==2
                    iAlias = find(~cellfun(@isempty,regexpi(scanID,imPar.tokenScanAliases(:,1),'once')));
                    if ~isempty(iAlias)
                        scanNames{iScan} = imPar.tokenScanAliases{iAlias,2};
                    else
                        % keep the original name
                        fprintf('No matching scan aliases found, keeping the original name...\n');
                        WarningMessage = ['ExploreASL_Import: Unknown scan ID ' scanID ' found, don"t know what this is'];
                        scan_name = scanNames{iScan};
                        branch = matches{1}; % Fallback (possibly not correct)
                        scanpath = fullfile(imPar.RawRoot,branch);
                        destdir = fullfile(SubjDir, [scanNames{iScan} '_' num2str(iSession)]); % Fallback (possibly not correct)
                        dcm2niiCatchedErrors = xASL_imp_CatchErrors('isempty(iAlias)', WarningMessage, dbstack, mfilename, pwd, scan_name, scanpath, destdir, dcm2niiCatchedErrors, imPar);
                    end
                end
                scan_name = scanNames{iScan};

                %% 4.3 Minimalistic feedback of where we are
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

                %% 4.4 Now pick the matching one from the folder list
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

                %% 4.5 Determine input and output paths
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
                
                %% 4.6 Start the conversion if this scan should not be skipped
                [imPar, globalCounts, x, summary_line, destdir, scanpath, scan_name, dcm2niiCatchedErrors, nii_files, first_match] = ...
                    xASL_imp_DCM2NII_Subject_StartConversion(...
                    imPar, globalCounts, x, bSkipThisOne, summary_line, destdir, scanpath, scan_name, dcm2niiCatchedErrors);
                
                
                %% 4.7 Store JSON files
                if ~isempty(nii_files)
                    jsonFiles = nii_files;
                    for iFile = 1:length(nii_files)
                        [Fpath, Ffile, ~] = fileparts(nii_files{iFile});
                        newFile = xASL_adm_GetFileList(Fpath,['^' Ffile '.json$'],'FPList');
                        jsonFiles{iFile} = newFile{1};
                    end
                    [parms, pathDcmDict] = xASL_imp_DCM2NII_Subject_StoreJSON(imPar, jsonFiles, first_match, settings.bUseDCMTK, pathDcmDict);
                end
                
                %% 4.8 In case of a single NII ASL file loaded from PAR/REC, we need to shuffle the dynamics from CCCC...LLLL order to CLCLCLCL... order
                [nii_files, summary_line, globalCounts, ASLContext] = xASL_imp_DCM2NII_Subject_ShuffleTheDynamics(globalCounts, scanpath, scan_name, nii_files, iSubject, iSession, iScan);
                
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
                
                %% 4.9 Make a copy of analysisdir in sourcedir
                xASL_imp_DCM2NII_Subject_CopyAnalysisDir(nii_files, settings.bClone2Source)
                
                % Copy single dicom as QC placeholder
                if settings.bCopySingleDicoms && ~isempty(first_match)
                    xASL_Copy(first_match, fullfile(destdir, ['DummyDicom_' scan_name '.dcm']), imPar.bOverwrite, imPar.bVerbose);
                end
                
                %% 4.10 Store the summary info so it can be sorted and printed below
                summary_lines{iSubject, iVisit, iSession, iScan} = summary_line;
            end % scansIDs
        end % sessionIDs
    end % visitIDs

        
end





