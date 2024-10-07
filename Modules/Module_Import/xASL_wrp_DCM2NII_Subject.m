function [x, PrintDICOMFields, dcm2niiCatchedErrors] = xASL_wrp_DCM2NII_Subject(x, matches, dcm2niiCatchedErrors)
%xASL_wrp_DCM2NII_Subject Run DCM2NII for one individual subject.
%
% FORMAT: [x, PrintDICOMFields, dcm2niiCatchedErrors] = xASL_wrp_DCM2NII_Subject(x, matches, dcm2niiCatchedErrors)
% 
% INPUT:
%   x                      - ExploreASL x structure (REQUIRED, STRUCT)
%   x.modules.import.imPar - Structure with import parameters (REQUIRED, STRUCT)
%   matches                - Matches (REQUIRED, CELL ARRAY)
%   dcm2niiCatchedErrors   - DCM2NII catched errors (REQUIRED, STRUCT)
%
% OUTPUT:
%   x                      - ExploreASL x structure (REQUIRED, STRUCT)
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
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     [x, PrintDICOMFields, dcm2niiCatchedErrors] = xASL_wrp_DCM2NII_Subject(x, matches, dcm2niiCatchedErrors);
% __________________________________
% Copyright 2015-2024 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________



    %% 1. Run DCM2NII for one individual subject
    
    % We do not iterate over subjects anymore, since this is done in xASL_Iteration now
    iSubject = find(strcmp(x.SUBJECT,x.SUBJECTS));
    
    % Overview subjects
    overviewSubjects = fieldnames(x.importOverview);
    
    % Get current subject
    thisSubject = x.importOverview.(overviewSubjects{iSubject});
    subjectID = x.modules.import.listsIDs.subjectIDs{iSubject};
    
    % Check subjectID
    thisSubject.subjectExport = xASL_imp_SubjectName(subjectID);
    
	x.modules.import.imPar.visitNames = thisSubject.visitIDs;

    %% 2. Iterate over visits
    for iVisit=1:thisSubject.nVisits
        
        % Get the current visit and the visit ID
        vFieldName = ['visit_' num2str(iVisit,'%03.f')];
        thisVisit = thisSubject.(vFieldName);
        
        % Determine the subject directory
        x.modules.import.SubjDir = xASL_imp_GetSubjDir(x, thisSubject.subjectExport, iVisit);

        if x.modules.import.imPar.SkipSubjectIfExists && exist(x.modules.import.SubjDir, 'dir')
            % we found the subject dir (i.e. SubjectVisit), so we skip it
            % this is ignored when x.modules.import.imPar.SkipSubjectIfExists is set to
            % false (default)
            continue
        end
        
        % Display subject-visit ID and add lock dir
		xASL_adm_BreakString('');
        fprintf('Importing subject = %s:   \n', [thisSubject.subjectExport '_' x.modules.import.imPar.visitNames{iVisit}]);

        %% 3. Loop over all sessions
        
        % Catch empty sessions (initial T1w e.g.)
        emptySessions = sum(cellfun(@isempty,thisVisit.runs))>0;
        indexEmptySession = -1; % Fallback
        if emptySessions
            indexEmptySession = find(cellfun(@isempty,thisVisit.runs));
        end
        
        for iSession=1:thisVisit.nSessions
            
            % Get current run
            vSessionName = ['run_' num2str(iSession,'%03.f')];
            thisRun = thisVisit.(vSessionName);

            % Get current run name
            x.modules.import.imPar.sessionNames{iSession} = thisRun.name;
            
            % Find empty sessions
            scanFields.emptySession = false;
            if sum(iSession==indexEmptySession)
                scanFields.emptySession = true;
            end
            
            % Quick & dirty fix for now...
            if ~(numel(thisRun.scanIDs)==numel(thisRun.ids))
                firstElement = {thisRun.ids{1}};
                thisRun.ids = cell(1,numel(thisRun.scanIDs));
                thisRun.ids(:) = firstElement;
            end
            
            %% 4. Iterate over scans
			% Safety-check, we can't allow multiple matches for the same sequence
			if numel(thisRun.scanIDs) ~= numel(unique(thisRun.scanIDs))
				fprintf('Detected duplicate scans: %s\n',strjoin(thisRun.scanIDs, ', '));
				error('We do not allow multiple token matches per scantype. Please adjust sourcestructure.json');
			end

            for iScan=1:numel(thisRun.scanIDs)
                % Get scan ID
                scanFields.runID = thisRun.ids{iScan};
                % Pack other scan related fields
                scanFields.subjectID = subjectID;
                scanFields.visitID = thisSubject.visitIDs{iVisit};
                scanFields.iSubject = iSubject;
                scanFields.iVisit = iVisit;
                scanFields.iSession = iSession;
                scanFields.iScan = iScan;
                scanFields.name = thisRun.name;

				% Create the session number based on the session name in format ASL_X
				if isempty(thisRun.name) || isempty(regexpi(thisRun.name, '^ASL_\d+$', 'once'))
					% Set to 1 if session name cannot be identified
 					scanFields.numberSession = 1;
				else
					scanFields.numberSession = xASL_str2num(thisRun.name(5:end));
				end
                % Convert scan
                [x, thisSubject,dcm2niiCatchedErrors, PrintDICOMFields] = ...
                    xASL_imp_DCM2NII_ConvertScan(x, matches, thisSubject, dcm2niiCatchedErrors, thisVisit, thisRun, scanFields);
            end
            
        end
        thisSubject.(vFieldName) = thisVisit;
        
    end
    
    % Make sure to update the subject id if there were illegal characters (according to BIDS)
    thisSubject.name = thisSubject.subjectExport;
    thisSubject = rmfield(thisSubject,'subjectExport');
    
    % Put data back into x structure
    x.importOverview.(overviewSubjects{iSubject}) = thisSubject;
    
        
end


%% Get subject directory
function SubjDir = xASL_imp_GetSubjDir(x, subjectExport, iVisit)

    % Only pad VisitID _1 _2 _3 etc if there are visits specified. Multiple visits is defined by the tokenVisitAliases.
    % If this is non-existing, it is set to 1, and if it does exist, it will put the _1 _2 _3 etc in the folder.
    % This fix allows to import a single visit from a range of specified visits.

    if x.modules.import.settings.bUseVisits
        if ~isempty(x.modules.import.imPar.visitNames{iVisit}) && strcmp(x.modules.import.imPar.visitNames{iVisit}(1),'_')
            % Only add '_' if there isn't one already
            SubjDir = fullfile(x.modules.import.imPar.TempRoot, [subjectExport x.modules.import.imPar.visitNames{iVisit}]);
        else
            % Subject/session directory with '_'
            SubjDir = fullfile(x.modules.import.imPar.TempRoot, [subjectExport '_' x.modules.import.imPar.visitNames{iVisit}]);
        end
    else
        SubjDir = fullfile(x.modules.import.imPar.TempRoot, subjectExport);
    end

end


%% If a subject ID has a `_` or a `-` or something else, we need to get rid of it
function subjectExport = xASL_imp_SubjectName(subjectID)

    % If a subject ID has a `_` or a `-` or something else, we need to get
    % rid of it in the directory name for BIDS, for the import of the
    % sourcedata we should still know about that though, which is why we
    % use a new variable for the subject ID directory export.
    
    subjectExport = xASL_adm_CorrectName(subjectID, 2);
    
    if ~strcmp(subjectID,subjectExport)
        fprintf(2,'Special characters in subject ID, changing %s to %s\n',subjectID,subjectExport);
    end
    
    % Since following import sub-modules depend on the x.importOverview field and
    % this is also based on the original read-only sourcedata, we have to
    % fix the subject ids in there afterwards, too!

end



