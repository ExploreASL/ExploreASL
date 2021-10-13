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


