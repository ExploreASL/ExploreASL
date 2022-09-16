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
if ~isfield(studyParFull,'ImportContexts')
	studyParSpecific = studyParFull;
	return;
end

% Initialize empty studyPar and then go through different contexts
studyParSpecific = struct();

% At each step, check for all conditions and set a flag for overwrite (first step having the lowest priority
for iContext = 1:length(studyParFull.ImportContexts)
	bOverwrite = 1;
	
	% If AliasHierarchy field is not there, then always take it
	if isfield(studyParFull.ImportContexts{iContext},'AliasHierarchy')
		aliasHierarchy = studyParFull.ImportContexts{iContext}.AliasHierarchy;
		
		% Check subject name - missing rule or missing regexp means that it is OK
		if ~isempty(subjectName) && ~isempty(aliasHierarchy) && ~isempty(aliasHierarchy{1})
			% If the regexp is not found, then we don't take this subejct
			if isempty(regexpi(subjectName,aliasHierarchy{1}))
				bOverwrite = 0;
			end
		end
		
		% Check sessions
		if ~isempty(sessionName) && length(aliasHierarchy)>1 && ~isempty(aliasHierarchy{2})
			if isempty(regexpi(sessionName,aliasHierarchy{2}))
				bOverwrite = 0;
			end
		end
		
		% Check runs
		if ~isempty(runName) && length(aliasHierarchy)>2 && ~isempty(aliasHierarchy{3})
			if isempty(regexpi(runName,aliasHierarchy{3}))
				bOverwrite = 0;
			end
		end
		
		% Remove the aliasHierarchy from this context
		studyParFull.ImportContexts{iContext} = rmfield(studyParFull.ImportContexts{iContext},'AliasHierarchy');	
	end
	
	if bOverwrite
		% Overwrite all fields
		listFields = fieldnames(studyParFull.ImportContexts{iContext});
		
		for iField = 1:length(listFields)
			studyParSpecific.(listFields{iField}) = studyParFull.ImportContexts{iContext}.(listFields{iField});
		end
	end
end
end
