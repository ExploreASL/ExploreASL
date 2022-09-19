function studyParSpecific = xASL_imp_StudyParPriority(studyParFull, subjectName, sessionName, runName)
%xASL_imp_StudyParPriority Takes the studyParFull a prioritizes based on the current subject/session/run
%
% FORMAT: studyParSpecific = xASL_imp_StudyParPriority(studyParFull[, subjectName, sessionName, runName])
%
% INPUT:
%   studyParFull - StudyPar possibly containing several contexts (REQUIRED, STRUCT)
%   subjectName  - Name of the current subject (OPTIONAL, DEFAULT = '')
%   sessionName  - Name of the current session (OPTIONAL, DEFAULT = '')
%   runName      - Name of the current run (OPTIONAL, DEFAULT = '')
%
% OUTPUT:
%   studyParSpecific - Resolved studyPar for a specific subject/session/run
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Takes studyPar with possibly several contexts and resolves the priority and saves the individual studyPar.
%                 First context has the lowest priority, following ones are more important and review previous variables. 
%                 The alias hierarchy within each context contains a list of regexp to be matched against subject/session/run.
%                 If subject doesn't exist or regexp doesn't exist or aliasHierarchy doesn't exist, it allows it.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        studyPar = xASL_imp_StudyParPriority(studyParFull, 'sub-01', '1', '2')
% __________________________________
% Copyright 2015-2022 ExploreASL

% Check existence of studyParFull and return if single context
if nargin < 1 || isempty(studyParFull)
	studyParSpecific = struct();
	return;
end

% Do nothing for a single-context studyPar
if ~isfield(studyParFull,'StudyPars')
	studyParSpecific = studyParFull;
	return;
end

% Initialize empty studyPar and then go through different contexts
studyParSpecific = struct();

% At each step, check for all conditions and set a flag for overwrite (first step having the lowest priority
for iContext = 1:length(studyParFull.StudyPars)
	bOverwrite = 1;
	
	% Checks for fields SubjectRegexp, VisitRegexp, SessionRegexp - if missing or empty, then take this context into account
	% Otherwise check if the respective field and Subject/Visit/Session name matches the regular expression
	
	% Check subject name - missing rule or missing regexp means that it is OK
	if isfield(studyParFull.StudyPars{iContext},'SubjectRegexp')
		if ~isempty(subjectName) && ~isempty(studyParFull.StudyPars{iContext}.SubjectRegexp)
			% If the regexp is not found, then we don't take this subejct
			if isempty(regexpi(subjectName, xASL_num2str(studyParFull.StudyPars{iContext}.SubjectRegexp)))
				bOverwrite = 0;
			end
		end
		studyParFull.StudyPars{iContext} = rmfield(studyParFull.StudyPars{iContext},'SubjectRegexp');	
	end
	
	% Check sessions
	if isfield(studyParFull.StudyPars{iContext},'VisitRegexp')
		if ~isempty(sessionName) && ~isempty(studyParFull.StudyPars{iContext}.VisitRegexp)
			% If the regexp is not found, then we don't take this subejct
			if isempty(regexpi(sessionName, xASL_num2str(studyParFull.StudyPars{iContext}.VisitRegexp)))
				bOverwrite = 0;
			end
		end
		studyParFull.StudyPars{iContext} = rmfield(studyParFull.StudyPars{iContext},'VisitRegexp');	
	end
		
	% Check runs
	if isfield(studyParFull.StudyPars{iContext},'SessionRegexp')
		if ~isempty(runName) && ~isempty(studyParFull.StudyPars{iContext}.SessionRegexp)
			% If the regexp is not found, then we don't take this subejct
			if isempty(regexpi(runName, xASL_num2str(studyParFull.StudyPars{iContext}.SessionRegexp)))
				bOverwrite = 0;
			end
		end
		studyParFull.StudyPars{iContext} = rmfield(studyParFull.StudyPars{iContext},'SessionRegexp');	
	end
	
	if bOverwrite
		% Overwrite all fields
		listFields = fieldnames(studyParFull.StudyPars{iContext});
		
		for iField = 1:length(listFields)
			studyParSpecific.(listFields{iField}) = studyParFull.StudyPars{iContext}.(listFields{iField});
		end
	end
end
end
