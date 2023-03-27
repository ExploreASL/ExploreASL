function [x, thisSubject, dcm2niiCatchedErrors, PrintDICOMFields] = xASL_imp_DCM2NII_ConvertScan(x, matches, thisSubject, dcm2niiCatchedErrors, thisVisit, thisRun, scanFields)
%xASL_imp_DCM2NII_ConvertScan Run DCM2NII for one individual scan.
%
% FORMAT: [x, thisSubject, dcm2niiCatchedErrors, PrintDICOMFields] = xASL_imp_DCM2NII_ConvertScan(x, matches, thisSubject, dcm2niiCatchedErrors, thisVisit, thisRun, scanFields)
% 
% INPUT:
%   x                      - ExploreASL x structure (REQUIRED, STRUCT)
%   x.modules.import.imPar - Structure with import parameters (REQUIRED, STRUCT)
%
% OUTPUT:
%   x                      - ExploreASL x structure (REQUIRED, STRUCT)
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Run DCM2NII for one individual scan.
%
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
% EXAMPLE:     n/a
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% 4.1 Initialize variables
    
    PrintDICOMFields = [];
    
    % Get other basic fields
    scanID = thisRun.scanIDs{scanFields.iScan};
    summary_line = [];
    thisSubject.summary_lines{scanFields.iSubject,scanFields.iVisit,scanFields.iSession,scanFields.iScan} = 'n/a';
    
    % If not verbose, track percentage progress
    if ~x.modules.import.imPar.bVerbose
        CounterN = (scanFields.iSession-1)*thisVisit.nScans+scanFields.iScan;
        CounterT = thisVisit.nSessions*thisVisit.nScans;
        xASL_TrackProgress(CounterN, CounterT);
    end

    %% 4.2 Convert scan ID to a suitable name and set scan-specific parameters
    if size(x.modules.import.imPar.tokenScanAliases,2)==2
        iAlias = find(~cellfun(@isempty,regexpi(scanID,x.modules.import.imPar.tokenScanAliases(:,1),'once')));
        if ~isempty(iAlias)
            thisVisit.scanNames{scanFields.iScan} = x.modules.import.imPar.tokenScanAliases{iAlias,2};
        else
            % keep the original name
            fprintf('No matching scan aliases found, keeping the original name...\n');
            WarningMessage = ['ExploreASL_Import: Unknown scan ID ' scanID ' found, don"t know what this is'];
            scan_name = thisVisit.scanNames{scanFields.iScan};
            % Fallback to first match (could be incorrect)
            branch = matches{1};
            scanpath = fullfile(x.modules.import.imPar.RawRoot,branch);
            % Determine destination directory
            destdir = fullfile(x.modules.import.SubjDir, [thisVisit.scanNames{scanFields.iScan} '_' num2str(scanFields.iSession)]);
            % Add fields to catched errors
            dcm2niiCatchedErrors = xASL_imp_CatchErrors('isempty(iAlias)', WarningMessage, dbstack, ...
                mfilename, pwd, scan_name, scanpath, destdir, dcm2niiCatchedErrors, x.modules.import.imPar);
        end
    end
    scan_name = thisVisit.scanNames{scanFields.iScan};

    %% 4.3 Minimalistic feedback of where we are
    if x.modules.import.imPar.bVerbose
        if isempty(scanFields.name)
            printSession = 'empty';
        else
            printSession = scanFields.name;
        end
        fprintf('Subject = %s, visit = %s, session = %s, scan = %s\n',scanFields.subjectID, scanFields.visitID, printSession, scan_name);
    end

%     % Defaults
    bOneScanIsEnough = false; % PM: THIS USED TO BE TRUE FOR ANATOMICAL SCANS
    bPutInSessionFolder = true; % THIS USED TO BE FALSE FOR ANATOMICAL SCANS
    if ~isempty(regexpi(scan_name, '^ASL4D.*')) || ~isempty(regexpi(scan_name, '^M0.*')) || ~isempty(regexpi(scan_name, '^func_bold$'))
		bPutInSessionFolder = true;
		bAnatomical = false;
	elseif ~isempty(regexpi(scan_name, '^T1.+')) || ~isempty(regexpi(scan_name, '^T2.+')) || ~isempty(regexpi(scan_name, '^FLAIR$')) || ~isempty(regexpi(scan_name, '^WMH_SEGM$'))
		%             bPutInSessionFolder = false;
		bAnatomical = true;
	else
		error('Unrecognized scantype');
    end

    
    %% 4.4 Now pick the matching one from the folder list
    
    % Only get the matching session
    iMatch = find(strcmp(x.modules.import.listsIDs.vSubjectIDs,scanFields.subjectID) & ...
        strcmp(x.modules.import.listsIDs.vVisitIDs, xASL_adm_CorrectName(scanFields.visitID,2,'_')) & ...
        strcmp(x.modules.import.listsIDs.vSessionIDs,scanFields.runID) & ...
        strcmpi(x.modules.import.listsIDs.vScanIDs,scanID));
    
    % Only report as missing if we need a scan for each session (i.e. ASL)
    if isempty(iMatch)
        if sum(thisSubject.globalCounts.converted_scans(scanFields.iSubject,scanFields.iVisit,:,scanFields.iScan))==0
            % Define warning message
            WarningMessage = ['Missing scan: ' [scanFields.subjectID x.modules.import.imPar.visitNames{scanFields.iVisit}] ', ' ...
                num2str(scanFields.iSession) ', ' scan_name];
            % Print warning
            if x.modules.import.imPar.bVerbose
                warning(WarningMessage);
            end
            % Add missing scans to global counts struct
            thisSubject.globalCounts.missing_scans(scanFields.iSubject, scanFields.iVisit, scanFields.iSession, scanFields.iScan) = 1;
        end
        % Add information to summary lines and stop import of this scan here
        thisSubject.summary_lines{scanFields.iSubject, scanFields.iVisit, scanFields.iSession, scanFields.iScan} = summary_line;
        return
    end

    %% 4.5 Determine input and output paths
    bSkipThisOne = false;
    branch = matches{iMatch};
    scanpath = fullfile(x.modules.import.imPar.RawRoot,branch);

    if ~isempty(strfind(thisVisit.scanNames{scanFields.iScan}, 'ASL4D')) || ~isempty(strfind(thisVisit.scanNames{scanFields.iScan}, 'M0'))
        session_name = scanFields.name;
    elseif ~isempty(strfind(thisVisit.scanNames{scanFields.iScan}, 'DSC4D'))
        session_name = ['DSC_' num2str(scanFields.iSession)];
    else
        % Allow multiple ScanTypes for sessions
        session_name = [thisVisit.scanNames{scanFields.iScan} '_' num2str(scanFields.iSession)]; 
    end

    % Determine the subject directory name in the temp folder
    if bPutInSessionFolder
        % Put in a subject_session folder
        destdir = fullfile(x.modules.import.SubjDir, session_name);
    else
        % Put in subject folder instead of session folder
        destdir = x.modules.import.SubjDir;
    end
    
    % One scan is enough, so skip this one if there was already a scan converted of this type (i.e. T1)
    if bOneScanIsEnough && sum(thisSubject.globalCounts.converted_scans(scanFields.iSubject,scanFields.iVisit,:,scanFields.iScan))~=0
        if x.modules.import.imPar.bVerbose
            fprintf('Skipping scan: %s, %s, %s\n',[scanFields.subjectID, x.modules.import.imPar.visitNames{scanFields.iVisit}], session_name, scan_name);
        end
        bSkipThisOne = true;
        destdir = [];
    end

    %% 4.6 Start the conversion if this scan should not be skipped
    
    % In case a previous run crashed we should remove the existing temp data here
    % (this does not really work ATM, since T1 files can be in the top directory of ASL e.g.)
    % if xASL_exist(destdir,'dir')
    %    fprintf(2,'Remove existing temp data...\n');
    %    xASL_delete(destdir,true);
    % end

    % Conversion
    [thisSubject.globalCounts, x, ~, destdir, scanpath, scan_name, dcm2niiCatchedErrors, nii_files, first_match] = ...
        xASL_imp_DCM2NII_Subject_StartConversion(...
        thisSubject.globalCounts, x, bSkipThisOne, summary_line, destdir, scanpath, scan_name, dcm2niiCatchedErrors);


    %% 4.7 Store JSON files
    jsonFiles = {};
    if ~isempty(nii_files)

        for iFile = 1:length(nii_files)
            [Fpath, Ffile, ~] = xASL_fileparts(nii_files{iFile});
            newFile = xASL_adm_GetFileList(Fpath,['^' Ffile '.json$'],'FPList');
            if ~isempty(newFile)
                jsonFiles{iFile} = newFile{1};
            else
                warning(['JSON sidecar missing for ' nii_files{iFile}]);
            end
        end
        if ~isempty(jsonFiles) && isempty(regexp(first_match,'[nii|nii\.gz]$'))
            [parms, x.modules.import.pathDcmDict] = xASL_imp_DCM2NII_Subject_StoreJSON(...
                x.modules.import.imPar, jsonFiles, first_match, x.modules.import.settings.bUseDCMTK, x.modules.import.pathDcmDict);
        end
    end

    
    %% 4.8 Sort ASL Volumes
    if ~bAnatomical % skip this for anatomical scans
        % In case of GE or PARREC/a single NII ASL file loaded from PAR/REC, we need to 
        % shuffle the dynamics from CCCC...LLLL order to CLCLCLCL... order.
        [x,nii_files, summary_line, thisSubject.globalCounts] = ...
            xASL_imp_DCM2NII_Subject_SortASLVolumes(x, thisSubject.globalCounts, scanpath, scan_name, nii_files, scanFields.iSubject, scanFields.iVisit, scanFields.iSession, scanFields.iScan);

        % Correct nifti rescale slope if parms.RescaleSlopeOriginal =~1 but nii.dat.scl_slope==1 (this can happen in case
        % of hidden scale slopes in private Philips header, that is dealt with by xASL_bids_Dicom2JSON but not by dcm2niiX.
        if ~isempty(nii_files) && exist('parms','var')
            [TempLine, PrintDICOMFields] = xASL_imp_AppendParmsParameters(parms);
            summary_line = [summary_line TempLine];
        end
    end

    %% 4.9 Copy single dicom as QC placeholder
    if x.modules.import.settings.bCopySingleDicoms && ~isempty(first_match) && ~regexpi(fExt, '\.(nii|nii\.gz)')
        xASL_Copy(first_match, fullfile(destdir, ['DummyDicom_' scan_name '.dcm']), x.modules.import.imPar.bOverwrite, x.modules.import.imPar.bVerbose);
    end

    %% 4.10 Store the summary info so it can be sorted and printed below
    thisSubject.summary_lines{scanFields.iSubject, scanFields.iVisit, scanFields.iSession, scanFields.iScan} = summary_line;

    %% 4.11 Check DCM2NIIX output
    xASL_imp_Check_DCM2NII_Output(nii_files,scanID);


end




%% Check the DCM2NII output
function xASL_imp_Check_DCM2NII_Output(nii_files,scanID)

    % For some ADNI cases there are multiple anatomical scans in a single session (case 006_S_4485 e.g.).
    % This can be troublesome for NII2BIDS and BIDS2LEGACY.
    if (~isempty(regexpi(scanID,'t1w')) ||  ~isempty(regexpi(scanID,'flair'))) && numel(nii_files)>1
        fprintf('Multiple anatomical NIfTIs for a single session...\n');
    end

end


