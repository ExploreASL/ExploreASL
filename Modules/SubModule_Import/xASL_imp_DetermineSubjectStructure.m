function x = xASL_imp_DetermineSubjectStructure(x)
%xASL_imp_DetermineSubjectStructure Determine subject/session/run structure from sourcedata or temp data
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
% DESCRIPTION:    Determine subject/session/run structure from sourcedata or temp data.
%         
%          The main goal is to create a sub-structure of x called x.overview.
%          x.overview does include a separate field for each subject.
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
% Copyright 2015-2021 ExploreASL


    %% Shared imPar related initialization
    if x.opts.ImportModules(1) || x.opts.ImportModules(2) || x.opts.ImportModules(3) || x.opts.ImportModules(4)
        if isfield(x.modules.import,'imPar') && isstruct(x.modules.import.imPar)
            % Basic import checks before execution
            x = xASL_imp_CheckImportSettings(x);
            % Check directories and permissions of sourcedata
            x = xASL_imp_CheckDirectoriesAndPermissions(x);
        else
            error('The imPar struct does not exist...');
        end
    else
        % Not a valid option
        error('Invalid option of import settings...');
    end


    %% Specific initialization for sourcedata, temp data, and rawdata
    if x.opts.ImportModules(1) && ...
            isfield(x.modules.import,'imPar') && isstruct(x.modules.import.imPar)
        % Determine structure from sourcedata
        x = xASL_imp_DetermineStructureFromSourcedata(x);
        
    elseif x.opts.ImportModules(2) && ...
            isfield(x.modules.import,'imPar') && isstruct(x.modules.import.imPar)
        % Determine structure from temp data
        x = xASL_imp_DetermineStructureFromTempdata(x);
        
    elseif (x.opts.ImportModules(3) || x.opts.ImportModules(4)) && ...
            isfield(x.modules.import,'imPar') && isstruct(x.modules.import.imPar)
        % Determine structure from rawdata
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
% Copyright 2015-2021 ExploreASL


    %% SUBJECTS
    x.SUBJECTS = xASL_adm_GetFileList(x.modules.import.imPar.BidsRoot,[],false,[],true);
    
    % Remove 'sub-' from subject name if it exists
    for iSubject=1:numel(x.SUBJECTS)
        if regexpi(x.SUBJECTS{iSubject},'sub-')==1
            x.SUBJECTS{iSubject} = x.SUBJECTS{iSubject}(length('sub-')+1:end);
        end 
    end
    
    if isempty(x.SUBJECTS)
        warning('Unable to find subjects in BIDS rawdata directory...');
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
% Copyright 2015-2021 ExploreASL


    %% Determine structure from temp data

    % Get subject/session list
    listSubjectsSessions = xASL_adm_GetFileList(x.modules.import.imPar.TempRoot,[],false,[],true);
    
    % Get subjects from list
    if numel(listSubjectsSessions)>0
        for iSubSes=1:numel(listSubjectsSessions)
            % Get current subject/session name
            curSubSes = listSubjectsSessions{iSubSes};
            % Determine subject name
            if ~isempty(regexp(curSubSes,'_','all')) && ~(numel(regexp(curSubSes,'_','all'))>1)
                % Multi-session notation
                x.SUBJECTS{iSubSes} = curSubSes(1:regexp(curSubSes,'_')-1);
            elseif isempty(regexp(curSubSes,'_','all'))
                % Single-session exceptions
                x.SUBJECTS{iSubSes} = curSubSes;
            else
                warning('It was not possible to determine the subject name from the temp data...');
                if numel(regexp(curSubSes,'_','all'))>1
                    warning('Multiple underscores in subject/session name...');
                end
            end
        end
    end
    
    % Unique the list
    x.SUBJECTS = unique(x.SUBJECTS);
    
    % Check if list is  empty
    if isempty(x.SUBJECTS)
        warning('Unable to find subjects in temp directory...');
    end
    


end


