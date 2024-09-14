function studyParSpecificSubjVisitSess = xASL_imp_StudyParPriority(studyParAll, subjectName, sessionName, runName, bVerbose)
%xASL_imp_StudyParPriority Takes the studyParAll a prioritizes based on the current subject/session/run
%
% FORMAT: studyParSpecificSubjVisitSess = xASL_imp_StudyParPriority(studyParAll[, subjectName, sessionName, runName, bVerbose])
%
% INPUT:
%   studyParAll  - StudyPar possibly containing several studyPar instances (REQUIRED, STRUCT)
%   subjectName  - Name of the current subject (OPTIONAL, DEFAULT = '')
%   sessionName  - Name of the current session (OPTIONAL, DEFAULT = '')
%   runName      - Name of the current run (OPTIONAL, DEFAULT = '')
%   bVerbose     - Write warnings when overwriting parameters (OPTIONAL, DEFAULT = true)
%
% OUTPUT:
%   studyParSpecificSubjVisitSess - Resolved studyPar for a specific subject/session/run
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
	studyParSpecificSubjVisitSess = struct();
	return;
end

% Default subject/visit/session strings to empty
if nargin < 2
	subjectName = '';
end

if nargin < 3
	sessionName = '';
end

if nargin < 4
	runName = '';
end

if nargin < 5 || isempty(bVerbose)
	bVerbose = true;
end

% Do nothing for a single-instance studyPar
if ~isfield(studyParAll,'StudyPars')
	studyParSpecificSubjVisitSess = studyParAll;
	return;
end

% Initialize empty studyPar and then go through different studyPars
studyParSpecificSubjVisitSess = struct();

textOverwriteWarning = '';

% At each step, check for all conditions and set a flag for overwrite (first step having the lowest priority
for iStudyPar = 1:length(studyParAll.StudyPars)
	bOverwrite = 1;
	
	% Checks for fields SubjectRegExp, VisitRegExp, SessionRegExp - if missing or empty, then take this specific studyPar into account
	% Otherwise check if the respective field and Subject/Visit/Session name matches the regular expression
	
	% Check subject name - missing rule or missing regexp means that the condition is fulfilled
	if isfield(studyParAll.StudyPars{iStudyPar},'SubjectRegExp')
		if ~isempty(subjectName) && ~isempty(studyParAll.StudyPars{iStudyPar}.SubjectRegExp)
			% If the regexp is not found, then we don't take this specific subject/visit/session
			if isempty(regexpi(subjectName, xASL_num2str(studyParAll.StudyPars{iStudyPar}.SubjectRegExp)))
				bOverwrite = 0;
			end
		end
	end
	
	% Check sessions
	if isfield(studyParAll.StudyPars{iStudyPar},'VisitRegExp')
		if ~isempty(sessionName) && ~isempty(studyParAll.StudyPars{iStudyPar}.VisitRegExp)
			% If the regexp is not found, then we don't take this specific subject/visit/session
			if isempty(regexpi(sessionName, xASL_num2str(studyParAll.StudyPars{iStudyPar}.VisitRegExp)))
				bOverwrite = 0;
			end
		end
	end
		
	% Check runs
	if isfield(studyParAll.StudyPars{iStudyPar},'SessionRegExp')
		if ~isempty(runName) && ~isempty(studyParAll.StudyPars{iStudyPar}.SessionRegExp)
			% If the regexp is not found, then we don't take this specific subject/visit/session 
			if isempty(regexpi(runName, xASL_num2str(studyParAll.StudyPars{iStudyPar}.SessionRegExp)))
				bOverwrite = 0;
			end
		end
	end
	
	if bOverwrite
		% Overwrite all fields
		listFields = fieldnames(studyParAll.StudyPars{iStudyPar});
		
		textOverwritingFields = '';
		for iField = 1:length(listFields)
			% list all fields that were overwritten
			if isfield(studyParSpecificSubjVisitSess, listFields{iField})
				textOverwritingFields = [textOverwritingFields ' ' listFields{iField} ','];
			end
			studyParSpecificSubjVisitSess.(listFields{iField}) = studyParAll.StudyPars{iStudyPar}.(listFields{iField});
		end
		
		% If there were overwritten fields, then we add the list of the fields and the RegExp list to the warning
		if ~isempty(textOverwritingFields)
			textOverwriteWarning = sprintf('%s For RegExp subject/visit/session ', textOverwriteWarning);
			if isfield(studyParAll.StudyPars{iStudyPar},'SubjectRegExp')
				textOverwriteWarning = [textOverwriteWarning xASL_num2str(studyParAll.StudyPars{iStudyPar}.SubjectRegExp) '/'];
			else
				textOverwriteWarning = [textOverwriteWarning  '/'];
			end
			if isfield(studyParAll.StudyPars{iStudyPar},'VisitRegExp')
				textOverwriteWarning = [textOverwriteWarning xASL_num2str(studyParAll.StudyPars{iStudyPar}.VisitRegExp) '/'];
			else
				textOverwriteWarning = [textOverwriteWarning  '/'];
			end
			if isfield(studyParAll.StudyPars{iStudyPar},'SessionRegExp')
				textOverwriteWarning = [textOverwriteWarning xASL_num2str(studyParAll.StudyPars{iStudyPar}.SessionRegExp)];
			end
			textOverwriteWarning = sprintf('%s, overwriting fields: %s\n', textOverwriteWarning, textOverwritingFields);
		end
	end
end

% If there are overwritten fields, then add a general info about the subject and print the warning
if ~isempty(textOverwriteWarning)
	textOverwriteWarning = sprintf('Evaluating multi-parameter studyPar for subject/visit/session: %s/%s/%s\n%s', subjectName, sessionName, runName, textOverwriteWarning);
	fprintf(textOverwriteWarning);
end

% Delete the keywords for merging
if isfield(studyParSpecificSubjVisitSess,'SubjectRegExp')
	studyParSpecificSubjVisitSess = rmfield(studyParSpecificSubjVisitSess,'SubjectRegExp');
end
if isfield(studyParSpecificSubjVisitSess,'VisitRegExp')
	studyParSpecificSubjVisitSess = rmfield(studyParSpecificSubjVisitSess,'VisitRegExp');
end
if isfield(studyParSpecificSubjVisitSess,'SessionRegExp')
	studyParSpecificSubjVisitSess = rmfield(studyParSpecificSubjVisitSess,'SessionRegExp');
end
end
