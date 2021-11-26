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
    
    % Check subjectID
    thisSubject.subjectExport = xASL_imp_SubjectName(subjectID);
    
    %% 2. Iterate over visits
    for iVisit=1:thisSubject.nVisits
        
        % Get the current visit and the visit ID
        vFieldName = ['visit_' num2str(iVisit,'%03.f')];
        thisVisit = thisSubject.(vFieldName);
        visitID = thisSubject.visitIDs{iVisit};

        % Convert visit ID to a suitable name
        if size(imPar.tokenVisitAliases,2)==2
            iAlias = find(~cellfun(@isempty,regexp(thisSubject.visitIDs{iVisit},imPar.tokenVisitAliases(:,1),'once')));
            if ~isempty(iAlias)
                imPar.visitNames{iVisit} = imPar.tokenVisitAliases{iAlias,2};
            end
        end
        
        % Determine the subject directory
        x.modules.import.SubjDir = xASL_imp_GetSubjDir(x,imPar,thisSubject.subjectExport,iVisit);

        if imPar.SkipSubjectIfExists && exist(x.modules.import.SubjDir, 'dir')
            % we found the subject dir (i.e. SubjectVisit), so we skip it
            % this is ignored when imPar.SkipSubjectIfExists is set to
            % false (default)
            continue
        end

        % Pad missing '_' if needed
        if ~strcmp(imPar.visitNames{iVisit}(1), '_')
            imPar.visitNames{iVisit} = ['_' imPar.visitNames{iVisit}];
        end
        
        % Display subject-visit ID and add lock dir
		xASL_adm_BreakString('');
        fprintf('Importing subject = %s:   \n', [thisSubject.subjectExport imPar.visitNames{iVisit}]);

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
            imPar.sessionNames{iSession} = thisRun.name;
            
            % Find empty sessions
            scanFields.emptySession = false;
            if iSession==indexEmptySession
                scanFields.emptySession = true;
            end
            
            % Quick & dirty fix for now...
            if ~(numel(thisRun.scanIDs)==numel(thisRun.ids))
                firstElement = {thisRun.ids{1}};
                thisRun.ids = cell(1,numel(thisRun.scanIDs));
                thisRun.ids(:) = firstElement;
            end
            
            %% 4. Iterate over scans
            for iScan=1:numel(thisRun.scanIDs)
                % Get scan ID
                scanFields.runID = thisRun.ids{iScan};
                % Pack other scan related fields
                scanFields.subjectID = subjectID;
                scanFields.visitID = visitID;
                scanFields.iSubject = iSubject;
                scanFields.iVisit = iVisit;
                scanFields.iSession = iSession;
                scanFields.iScan = iScan;
                scanFields.name = thisRun.name;
                % Convert scan
                [x,imPar,thisSubject,dcm2niiCatchedErrors,PrintDICOMFields] = ...
                    xASL_imp_DCM2NII_ConvertScan(x,imPar,matches,thisSubject,dcm2niiCatchedErrors,thisVisit,thisRun,scanFields);
            end
            
        end
        thisSubject.(vFieldName) = thisVisit;
        
    end
    
    % Make sure to update the subject id if there were illegal characters (according to BIDS)
    thisSubject.name = thisSubject.subjectExport;
    thisSubject = rmfield(thisSubject,'subjectExport');
    
    % Put data back into x structure
    x.overview.(overviewSubjects{iSubject}) = thisSubject;
    
        
end


%% Get subject directory
function SubjDir = xASL_imp_GetSubjDir(x, imPar, subjectExport, iVisit)

    % Only pad VisitID _1 _2 _3 etc if there are visits specified. Multiple visits is defined by the tokenVisitAliases.
    % If this is non-existing, it is set to 1, and if it does exist, it will put the _1 _2 _3 etc in the folder.
    % This fix allows to import a single visit from a range of specified visits.

    if x.modules.import.settings.bUseVisits
        if ~isempty(imPar.visitNames{iVisit}) && strcmp(imPar.visitNames{iVisit}(1),'_')
            % Only add '_' if there isn't one already
            SubjDir = fullfile(imPar.TempRoot, [subjectExport imPar.visitNames{iVisit}]);
        else
            % Subject/session directory with '_'
            SubjDir = fullfile(imPar.TempRoot, [subjectExport '_' imPar.visitNames{iVisit}]);
        end
    else
        SubjDir = fullfile(imPar.TempRoot, subjectExport);
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
    
    % Since following import sub-modules depend on the x.overview field and
    % this is also based on the original read-only sourcedata, we have to
    % fix the subject ids in there afterwards, too!

end



