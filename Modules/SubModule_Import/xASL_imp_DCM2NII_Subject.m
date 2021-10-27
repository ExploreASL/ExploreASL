function [x, imPar, PrintDICOMFields, dcm2niiCatchedErrors] = xASL_imp_DCM2NII_Subject(x, imPar, matches, dcm2niiCatchedErrors)
%xASL_imp_DCM2NII_Subject Run DCM2NII for one individual subject.
%
% FORMAT: [x, imPar, PrintDICOMFields, dcm2niiCatchedErrors] = xASL_imp_DCM2NII_Subject(x, imPar, matches, dcm2niiCatchedErrors)
% 
% INPUT:
%   x                      - ExploreASL x structure (REQUIRED, STRUCT)
%   imPar                  - Structure with import parameters (REQUIRED, STRUCT)
%   matches                - Matches (REQUIRED, CELL ARRAY)
%   dcm2niiCatchedErrors   - DCM2NII catched errors (REQUIRED, STRUCT)
%
% OUTPUT:
%   x                      - ExploreASL x structure (REQUIRED, STRUCT)
%   imPar                  - Structure with import parameters 
%   PrintDICOMFields       - Print DICOM fields
%   dcm2niiCatchedErrors   - DCM2NII catched errors
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Run DCM2NII for one individual subject.
%
% 1. Run DCM2NII for one individual subject
% 2. Iterate over visits
% 3. Loop through all sessions
% 4. Iterate over scans
% -  1. Initialize variables (scanID, summary_line, first_match)
% -  2. Convert scan ID to a suitable name and set scan-specific parameters
% -  3. Minimalistic feedback of where we are
% -  4. Now pick the matching one from the folder list
% -  5. Determine input and output paths
% -  6. Start the conversion if this scan should not be skipped
% -  7. Store JSON files
% -  8. In case of a single NII ASL file loaded from PAR/REC, we need to shuffle the dynamics from CCCC...LLLL order to CLCLCLCL... order
% -  9  Copy single dicom as QC placeholder
% - 10. Store the summary info so it can be sorted and printed below
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     [x, imPar, PrintDICOMFields, dcm2niiCatchedErrors] = xASL_imp_DCM2NII_Subject(x, imPar, matches, dcm2niiCatchedErrors);
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% 1. Run DCM2NII for one individual subject
    
    % We do not iterate over subjects anymore, since this is done in xASL_Iteration now
    iSubject = find(strcmp(x.SUBJECT,x.SUBJECTS));
    
    % Overview subjects
    overviewSubjects = fieldnames(x.overview);
    
    % Get current subject
    thisSubject = x.overview.(overviewSubjects{iSubject});
    subjectID = x.modules.import.listsIDs.subjectIDs{iSubject};
    
    %% 2. Iterate over visits
    for iVisit=1:thisSubject.nVisits
        
        % Get fieldname
        vFieldName = ['visit_' num2str(iVisit,'%03.f')];
        
        % Get visit
        thisVisit = thisSubject.(vFieldName);
        
        % Get visit ID
        visitID = thisSubject.visitIDs{iVisit};

        % Convert visit ID to a suitable name
        if size(imPar.tokenVisitAliases,2)==2
            iAlias = find(~cellfun(@isempty,regexp(thisSubject.visitIDs{iVisit},imPar.tokenVisitAliases(:,1),'once')));
            if ~isempty(iAlias)
                imPar.visitNames{iVisit} = imPar.tokenVisitAliases{iAlias,2};
            end
        end

        if x.modules.import.settings.bUseVisits 
            
            % Only pad VisitID _1 _2 _3 etc if there are visits specified
            % Multiple visits is defined by the tokenVisitAliases.
            % If this is non-existing, it is set to 1, and if it does exist, it will put the _1 _2 _3 etc in the folder
            % this fix allows to import a single visit from a range of specified visits.
            
            if ~isempty(imPar.visitNames{iVisit}) && strcmp(imPar.visitNames{iVisit}(1),'_')
                % Only add '_' if there isn't one already
                SubjDir = fullfile(imPar.TempRoot, [subjectID imPar.visitNames{iVisit}]);
            else
                % Subject/session directory with '_'
                SubjDir = fullfile(imPar.TempRoot, [subjectID '_' imPar.visitNames{iVisit}]);
            end
        else
            SubjDir = fullfile(imPar.TempRoot, subjectID);
        end

        if imPar.SkipSubjectIfExists && exist(SubjDir, 'dir')
            continue; % we found the subject dir (i.e. SubjectVisit), so we skip it
            % this is ignored when imPar.SkipSubjectIfExists is set to
            % false (default)
        end

        % Pad missing '_' if needed
        if ~strcmp(imPar.visitNames{iVisit}(1), '_')
            imPar.visitNames{iVisit} = ['_' imPar.visitNames{iVisit}];
        end
        
        % Display subject-visit ID and add lock dir
        fprintf('%s\nImporting subject=%s:   \n',...
            '==============================================================================================',...
            [subjectID imPar.visitNames{iVisit}]);

        %% 3. Loop through all sessions
        for iSession=1:thisVisit.nSessions
            
            % Get current run
            vSessionName = ['run_' num2str(iSession,'%03.f')];
            thisRun = thisVisit.(vSessionName);
            sessionID = thisRun.ids;

            % convert session ID to a suitable name
            imPar.sessionNames{iSession} = thisRun.name;
            
            %% 4. Iterate over scans
            for iScan=1:thisVisit.nScans
                
                %% 4.1 Initialize variables (scanID, summary_line, first_match)
                scanID = thisVisit.scanIDs{iScan};                         % I think here we actually need the scan IDs of the current run, not the current visit?!?
                summary_line = [];
                first_match = [];
                thisSubject.summary_lines{iSubject,iVisit,iSession,iScan} = 'n/a';

                if ~imPar.bVerbose % if not verbose, track % progress
                    CounterN = (iSession-1)*thisVisit.nScans+iScan;
                    CounterT = thisVisit.nSessions*thisVisit.nScans;
                    xASL_TrackProgress(CounterN, CounterT);
                end

                %% 4.2 Convert scan ID to a suitable name and set scan-specific parameters
                if size(imPar.tokenScanAliases,2)==2
                    iAlias = find(~cellfun(@isempty,regexpi(scanID,imPar.tokenScanAliases(:,1),'once')));
                    if ~isempty(iAlias)
                        thisVisit.scanNames{iScan} = imPar.tokenScanAliases{iAlias,2};
                    else
                        % keep the original name
                        fprintf('No matching scan aliases found, keeping the original name...\n');
                        WarningMessage = ['ExploreASL_Import: Unknown scan ID ' scanID ' found, don"t know what this is'];
                        scan_name = thisVisit.scanNames{iScan};
                        branch = matches{1}; % Fallback (possibly not correct)
                        scanpath = fullfile(imPar.RawRoot,branch);
                        destdir = fullfile(SubjDir, [thisVisit.scanNames{iScan} '_' num2str(iSession)]); % Fallback (possibly not correct)
                        dcm2niiCatchedErrors = xASL_imp_CatchErrors('isempty(iAlias)', WarningMessage, dbstack, mfilename, pwd, scan_name, scanpath, destdir, dcm2niiCatchedErrors, imPar);
                    end
                end
                scan_name = thisVisit.scanNames{iScan};

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
                iMatch = find(strcmp( x.modules.import.listsIDs.vSubjectIDs,subjectID) & ...
                    strcmp( x.modules.import.listsIDs.vVisitIDs, xASL_adm_CorrectName(visitID,2,'_')) & ...
                    strcmp( x.modules.import.listsIDs.vSessionIDs,sessionID) & ...
                    strcmpi( x.modules.import.listsIDs.vScanIDs,scanID) ); % only get the matching session
                if isempty(iMatch)
                    % only report as missing if we need a scan for each session (i.e. ASL)
                    if sum(thisSubject.globalCounts.converted_scans(iSubject,iVisit,:,iScan))==0
                        WarningMessage = ['Missing scan: ' [subjectID imPar.visitNames{iVisit}] ', ' num2str(iSession) ', ' scan_name];
                        if imPar.bVerbose; warning(WarningMessage); end
                        thisSubject.globalCounts.missing_scans(iSubject, iVisit, iSession, iScan) = 1;
                    end

                    thisSubject.summary_lines{iSubject, iVisit, iSession, iScan} = summary_line;
                    continue;
                    warning('Dont forget to comment continue here for debugging');
                end

                %% 4.5 Determine input and output paths
                bSkipThisOne = false;
                branch = matches{iMatch};
                scanpath = fullfile(imPar.RawRoot,branch);

                if ~isempty(strfind(thisVisit.scanNames{iScan}, 'ASL4D')) || ~isempty(strfind(thisVisit.scanNames{iScan}, 'M0'))
                    session_name = ['ASL_' num2str(iSession)];
                elseif ~isempty(strfind(thisVisit.scanNames{iScan}, 'DSC4D'))
                    session_name = ['DSC_' num2str(iSession)];
                else
                    session_name = [thisVisit.scanNames{iScan} '_' num2str(iSession)]; % Allow multiple ScanTypes for sessions
                end

                if bPutInSessionFolder
                    destdir = fullfile(SubjDir, session_name);
                else % put in subject folder instead of session folder
                    destdir = SubjDir;
                end

                if bOneScanIsEnough && sum(thisSubject.globalCounts.converted_scans(iSubject,iVisit,:,iScan))~=0
                    % one scan is enough, so skip this one if there was already a scan converted of this type (i.e. T1)
                    if imPar.bVerbose; fprintf('Skipping scan: %s, %s, %s\n',[subjectID imPar.visitNames{iVisit}],session_name,scan_name); end
                    bSkipThisOne = true;
                    destdir = []; % just in case
                end
                
                %% 4.6 Start the conversion if this scan should not be skipped
                [imPar, thisSubject.globalCounts, x, summary_line, destdir, scanpath, scan_name, dcm2niiCatchedErrors, nii_files, first_match] = ...
                    xASL_imp_DCM2NII_Subject_StartConversion(...
                    imPar, thisSubject.globalCounts, x, bSkipThisOne, summary_line, destdir, scanpath, scan_name, dcm2niiCatchedErrors);
                
                
                %% 4.7 Store JSON files
                if ~isempty(nii_files)
                    jsonFiles = nii_files;
                    for iFile = 1:length(nii_files)
                        [Fpath, Ffile, ~] = fileparts(nii_files{iFile});
                        newFile = xASL_adm_GetFileList(Fpath,['^' Ffile '.json$'],'FPList');
                        jsonFiles{iFile} = newFile{1};
                    end
                    [parms, x.modules.import.pathDcmDict] = xASL_imp_DCM2NII_Subject_StoreJSON(imPar, jsonFiles, first_match, x.modules.import.settings.bUseDCMTK, x.modules.import.pathDcmDict);
                end
                
                %% 4.8 In case of GE or PARREC/a single NII ASL file loaded from PAR/REC, we need to shuffle the dynamics from CCCC...LLLL order to CLCLCLCL... order
                [x,nii_files, summary_line, thisSubject.globalCounts] = ...
                    xASL_imp_DCM2NII_Subject_SortASLVolumes(x, thisSubject.globalCounts, scanpath, scan_name, nii_files, iSubject, iSession, iScan);

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
                
                %% 4.9 Copy single dicom as QC placeholder
                if x.modules.import.settings.bCopySingleDicoms && ~isempty(first_match)
                    xASL_Copy(first_match, fullfile(destdir, ['DummyDicom_' scan_name '.dcm']), imPar.bOverwrite, imPar.bVerbose);
                end
                
                %% 4.10 Store the summary info so it can be sorted and printed below
                thisSubject.summary_lines{iSubject, iVisit, iSession, iScan} = summary_line;
                
                %% 4.11 Check DCM2NIIX output
                xASL_imp_Check_DCM2NII_Output(nii_files,scanID);
                
                % scansIDs    
            end
            % sessionIDs
        end
        thisSubject.(vFieldName) = thisVisit;
        % visitIDs
    end
    
    
    %% Put data back into x structure
    x.overview.(overviewSubjects{iSubject}) = thisSubject;
        
end


%% Check the DCM2NII output
function xASL_imp_Check_DCM2NII_Output(nii_files,scanID)

    % For some ADNI cases there are multiple anatomical scans in a single session (case 006_S_4485 e.g.).
    % This can be troublesome for NII2BIDS and BIDS2LEGACY.
    if (~isempty(regexpi(scanID,'t1w')) ||  ~isempty(regexpi(scanID,'flair'))) && numel(nii_files)>1
        fprintf('Multiple anatomical NIfTIs for a single session...\n');
    end

end





