function x = xASL_imp_DetermineSubjectStructure(x)
%xASL_imp_DetermineSubjectStructure Determine subject/visit/session structure from sourcedata or temp data
%
% FORMAT: x = xASL_imp_DetermineSubjectStructure(x)
%
% INPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Determine subject/visit/session structure from sourcedata or temp data.
%         
%          The main goal is to create a sub-structure of x called x.importOverview.
%          x.importOverview does include a separate field for each subject.
%          Each subject has visit fields and visits have run fields.
%          All of this is used to track the number of subjects, number of visits, number of runs,
%          number of scans, etc. reliably. Within later parts of the pipeline it is also possible
%          to extract the current subject, current visit or current run to determine the part of
%          a population/subject/visit/run we're working on right now. We mostly used the variables
%          thisSubject, thisVisit, etc. for that.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright (c) 2015-2024 ExploreASL

    %% Check if imPar exists
    imParCondition = isfield(x.modules.import,'imPar') && isstruct(x.modules.import.imPar);

    %% Specific initialization for sourcedata, temp data, and rawdata
    if x.opts.bImport(1) && x.opts.bImportData && imParCondition
        % Determine structure from sourcedata
        x = xASL_imp_DetermineStructureFromSourcedata(x);
        
    elseif x.opts.bImport(2) && x.opts.bImportData && imParCondition
        % Determine structure from temp data
        x = xASL_imp_DetermineStructureFromTempdata(x);
        
    elseif x.opts.bImport(3)
        % Determine structure from rawdata
        warning('Loading subject-list from rawdata for defacing, this may differ from subjects in sourcedata');
        x = xASL_imp_DetermineStructureFromRawdata(x);
    end
    
    % SESSIONS DUMMY
    x.SESSIONS = {''};
end

function [x] = xASL_imp_DetermineStructureFromRawdata(x)
%xASL_imp_DetermineStructureFromRawdata Determine structure from rawdata
%
% FORMAT: [x] = xASL_imp_DetermineStructureFromRawdata(x)
%
% INPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Determine structure from rawdata.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright 2015-2024 ExploreASL

    %% Check if rawdata exists
    if xASL_exist(x.dir.RawData,'dir')
        
        % SUBJECTS
        x.SUBJECTS = xASL_adm_GetFileList(x.dir.RawData,[],false,[],true);

        % Remove 'sub-' from subject name if it exists
        for iSubject=1:numel(x.SUBJECTS)
            if regexpi(x.SUBJECTS{iSubject},'sub-')==1
                x.SUBJECTS{iSubject} = x.SUBJECTS{iSubject}(length('sub-')+1:end);
            end 
        end
   else
        % Maybe a user does not have a BIDS sourcedata or rawdata directory
        % and only runs the ExploreASL workflow on derivatives data.
        % Maybe BIDS2Legacy is turned on, but it actually shouldn't be.
        fprintf(2,'There is no rawdata directory, skipping defacing...\n');
    end
end

function [x] = xASL_imp_DetermineStructureFromTempdata(x)
%xASL_imp_DetermineStructureFromTempdata Determine structure from temp data
%
% FORMAT: [x] = xASL_imp_DetermineStructureFromTempdata(x)
%
% INPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Determine structure from temp data.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright 2015-2024 ExploreASL

    %% Determine structure from temp data

    % Get subject/visits list
    listSubjectsVisits = xASL_adm_GetFileList(x.modules.import.imPar.TempRoot,[],false,[],true);
    
	% Initialize an empty list of visits
	listVisits = {};
	
    % Get subjects from list
    if numel(listSubjectsVisits)>0
        for iSubVis=1:numel(listSubjectsVisits)
            % Get current subject/visit name
            curSubVis = listSubjectsVisits{iSubVis};
            % Determine subject name
			indexSeparator = regexp(curSubVis, '_', 'all');
            if isempty(indexSeparator)
                % Single-visit exceptions
                x.SUBJECTS{iSubVis} = curSubVis;
				% But we don't add a visit name to the list as no visit name was detected
			elseif numel(indexSeparator) == 1
                % Multi-visit notation
                x.SUBJECTS{iSubVis} = curSubVis(1:indexSeparator-1);
				% A visit name (non-empty but potentially also empty) was detected. We add it to the list of visits
				% The list of visits is only used here to construct tokenVisitAliases if missing, so the indexing of x.SUBJECTS and listVisits could differ
				listVisits{numel(listVisits)+1} = curSubVis(indexSeparator+1:end);
			else
				% Multiple separators
                warning('In BIDS, subject and session names shouold not contain hyphens (-) or underscores (_) as these are used as separators. Removing them in the final BIDS values/names');
            end
        end
    end
    
    % Unique the list
    if isfield(x, 'SUBJECTS')
        % Get the unique subjects
        x.SUBJECTS = unique(x.SUBJECTS);
        % Check if list is  empty
        if isempty(x.SUBJECTS)
            warning('Unable to find subjects in temp directory...');
        end
    else
        warning('x.SUBJECTS is undefined, data loading is going to fail...');
    end
    
	% If tokenVisitAliases is not provided, we assume visits are not renamed and create a dummy tokenVisitAliases list
	if ~isfield(x.modules.import.imPar, 'tokenVisitAliases') || isempty(x.modules.import.imPar.tokenVisitAliases)
		% Create a unique list of visits
		listVisits = unique(listVisits);

		% Create a list of unique visits names
		if isempty(listVisits)
			% In case we haven't detected any visits (not even with empty visit names), we set them to empty
			x.modules.import.imPar.tokenVisitAliases = [];
		else
			x.modules.import.imPar.tokenVisitAliases = listVisits(:);
			x.modules.import.imPar.tokenVisitAliases(:,2) = x.modules.import.imPar.tokenVisitAliases(:,1);
		end
	end
  
end
