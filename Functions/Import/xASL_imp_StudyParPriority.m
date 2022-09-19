function studyParSpecific = xASL_imp_StudyParPriority(studyParAll, subjectName, sessionName, runName)
%xASL_imp_StudyParPriority Takes the studyParAll a prioritizes based on the current subject/session/run
%
% FORMAT: studyParSpecific = xASL_imp_StudyParPriority(studyParAll[, subjectName, sessionName, runName])
%
% INPUT:
%   studyParAll  - StudyPar possibly containing several studyPar instances (REQUIRED, STRUCT)
%   subjectName  - Name of the current subject (OPTIONAL, DEFAULT = '')
%   sessionName  - Name of the current session (OPTIONAL, DEFAULT = '')
%   runName      - Name of the current run (OPTIONAL, DEFAULT = '')
%
% OUTPUT:
%   studyParSpecific - Resolved studyPar for a specific subject/session/run
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Takes studyPar with possibly several studyPar instances and resolves the priority and saves the individual studyPar.
%                 First studyPar instance has the lowest priority, following ones are more important and review previous variables. 
%                 The alias hierarchy within each studyPar instance contains a list of regexp to be matched against subject/session/run.
%                 If subject doesn't exist or regexp doesn't exist or aliasHierarchy doesn't exist, it allows it.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        studyPar = xASL_imp_StudyParPriority(studyParAll, 'sub-01', '1', '2')
% __________________________________
% Copyright 2015-2022 ExploreASL

% Check existence of studyParAll and return if single studyPar
if nargin < 1 || isempty(studyParAll)
	studyParSpecific = struct();
	return;
end

% Do nothing for a single-instance studyPar
if ~isfield(studyParAll,'StudyPars')
	studyParSpecific = studyParAll;
	return;
end

% Initialize empty studyPar and then go through different studyPars
studyParSpecific = struct();

% At each step, check for all conditions and set a flag for overwrite (first step having the lowest priority
for iStudyPar = 1:length(studyParAll.StudyPars)
	bOverwrite = 1;
	
	% Checks for fields SubjectRegexp, VisitRegexp, SessionRegexp - if missing or empty, then take this specific studyPar into account
	% Otherwise check if the respective field and Subject/Visit/Session name matches the regular expression
	
	% Check subject name - missing rule or missing regexp means that the condition is fulfilled
	if isfield(studyParAll.StudyPars{iStudyPar},'SubjectRegexp')
		if ~isempty(subjectName) && ~isempty(studyParAll.StudyPars{iStudyPar}.SubjectRegexp)
			% If the regexp is not found, then we don't take this specific subject/visit/session
			if isempty(regexpi(subjectName, xASL_num2str(studyParAll.StudyPars{iStudyPar}.SubjectRegexp)))
				bOverwrite = 0;
			end
		end
		studyParAll.StudyPars{iStudyPar} = rmfield(studyParAll.StudyPars{iStudyPar},'SubjectRegexp');	
	end
	
	% Check sessions
	if isfield(studyParAll.StudyPars{iStudyPar},'VisitRegexp')
		if ~isempty(sessionName) && ~isempty(studyParAll.StudyPars{iStudyPar}.VisitRegexp)
			% If the regexp is not found, then we don't take this specific subject/visit/session
			if isempty(regexpi(sessionName, xASL_num2str(studyParAll.StudyPars{iStudyPar}.VisitRegexp)))
				bOverwrite = 0;
			end
		end
		studyParAll.StudyPars{iStudyPar} = rmfield(studyParAll.StudyPars{iStudyPar},'VisitRegexp');	
	end
		
	% Check runs
	if isfield(studyParAll.StudyPars{iStudyPar},'SessionRegexp')
		if ~isempty(runName) && ~isempty(studyParAll.StudyPars{iStudyPar}.SessionRegexp)
			% If the regexp is not found, then we don't take this specific subject/visit/session 
			if isempty(regexpi(runName, xASL_num2str(studyParAll.StudyPars{iStudyPar}.SessionRegexp)))
				bOverwrite = 0;
			end
		end
		studyParAll.StudyPars{iStudyPar} = rmfield(studyParAll.StudyPars{iStudyPar},'SessionRegexp');	
	end
	
	if bOverwrite
		% Overwrite all fields
		listFields = fieldnames(studyParAll.StudyPars{iStudyPar});
		
		for iField = 1:length(listFields)
			studyParSpecific.(listFields{iField}) = studyParAll.StudyPars{iStudyPar}.(listFields{iField});
		end
	end
end
end
