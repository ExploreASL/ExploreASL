function x = xASL_init_PrintCheckSettings(x)
%xASL_init_PrintCheckSettings Check whether pre-defined settings existed in dataPar.json
%
% FORMAT: x = xASL_init_PrintCheckSettings(x)
%
% INPUT:
%   x       - ExploreASL x structure (STRUCT, REQUIRED)
%
% OUTPUT:
%   x       - ExploreASL x structure (STRUCT)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Check whether pre-defined settings existed in `dataPar.json`.
%
% Prints these on the screen as the start of the pipeline.
% Runs following steps:
%
% 1. Set default settings if not defined
% 2. Print data/study specific settings
% 3. Print warnings
%
% EXAMPLE:     This is part of the initialization workflow. Check out the usage there.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:  n/a
%
% Copyright 2015-2021 ExploreASL
    

    %% -----------------------------------------------------------------------
    %% 1) Set default settings if not defined
    if ~isfield(x.settings,'Quality') || (x.settings.Quality~=0 && x.settings.Quality~=1)
        x.settings.Quality = 1;
        fprintf('%s\n', 'Default Quality=1 used (optimal quality)');
    end
    if ~isfield(x.settings,'DELETETEMP') || (x.settings.DELETETEMP~=0 && x.settings.DELETETEMP~=1)
        x.settings.DELETETEMP = 1;
    %     fprintf('%s\n','Default x.settings.DELETETEMP=1 used (delete files temporarily used for processing)');
    end

    %% -----------------------------------------------------------------------
    %% 2) Print data/study specific settings
    xASL_adm_BreakString('Additional Settings',[],[],1);
    
    if x.opts.nWorkers>1
        fprintf(['I am worker ' num2str(x.opts.iWorker) '/' num2str(x.opts.nWorkers) '\n']);
        fprintf('Note that the resulting number of scans mentioned below applies only to this worker\n');
    end

    fprintf('%s\n',[num2str(x.dataset.nTotalSubjects) ' scans - ' ...
        num2str(x.dataset.nExcluded) ' exclusions, resulting in ' ...
        num2str(x.nSubjects) ' scans of: ']);

    for iT=1:x.dataset.nTimePointsTotal
        fprintf('%s\n',['Longitudinal timePoint ' num2str(iT) ' = ' ...
            num2str(x.dataset.nTimePointTotalSubjects(iT)) ' scans - ' ...
            num2str(x.dataset.nTimePointExcluded(iT)) ' exclusions = ' ...
            num2str(x.dataset.nTimePointSubjects(iT)) ' scans']);
    end

    fprintf('%s\n',['ASL sessions: ' num2str(x.dataset.nSessions)]);

    fprintf('\n%s','Ancillary data, sets: ');
    if isfield(x.S,'SetsID')
            fprintf('%s\n',[num2str(size(x.S.SetsID,2)) ' sets are defined for ' ...
                num2str(size(x.S.SetsID,1)) ' "SubjectsSessions"']);

            for iSet=1:size(x.S.SetsID,2)
                fprintf(['Set ' num2str(iSet) ' = "' x.S.SetsName{iSet} '" options ']);
                for iOption=1:size(x.S.SetsOptions{iSet},2)
                    fprintf(['"' x.S.SetsOptions{iSet}{iOption} '"']);
                    if iOption~= size(x.S.SetsOptions{iSet},2)
                        fprintf(' & ');
                    end
                end
                if      x.S.Sets1_2Sample(iSet)==1
                        fprintf(', codes for paired data');
                elseif  x.S.Sets1_2Sample(iSet)==2
                        fprintf(', codes for two-sample data');
                elseif  x.S.Sets1_2Sample(iSet)==3
                        fprintf([', continuous variate (with ' ...
                            num2str(length(unique(x.S.SetsID(:,iSet)))) ' unique values)']);
                end                    

                fprintf('\n');
            end
    else    
        fprintf('%s\n','No sets are defined');
    end

    if ~isfield(x.Q,'M0')
    %     warning('M0 option missing!');
    else
        fprintf('\n%s\n',['M0 option selected is "' num2str(x.Q.M0) '"']);
    end
    
    % Get the ExploreASL text width
    textBreakExploreASL = xASL_adm_BreakString([], [], false, 0, 0);

    if length(x.D.ROOT)>length(textBreakExploreASL)
        fprintf('x.D.ROOT  %s\n', x.D.ROOT(1:length(textBreakExploreASL)-10));
        fprintf('          ... %s\n', x.D.ROOT(length(textBreakExploreASL)-9:end));
    else
        fprintf('x.D.ROOT  %s\n', x.D.ROOT);
    end
    fprintf('x.settings.DELETETEMP %s\n',[num2str(x.settings.DELETETEMP) ' (delete temporary files)']);
    fprintf('x.settings.Quality    %s\n',[num2str(x.settings.Quality) ' (0 = fast try-out; 1 = normal high quality)']);


    %% -----------------------------------------------------------------------
    %% 3) Print warnings
    xASL_adm_BreakString('');
    fprintf('\n');
    field_symbol = {'subjectRegexp'};

    for iField=1:length(field_symbol)
        if ~isfield(x.dataset,field_symbol{iField})
            warning(['x.dataset' field_symbol{iField} ' was not defined in dataPar.json!'])
        end
    end

    if ~isfield(x,'D')
        warning('x.D didn''nt exist');
    else
        field_symbol = {'ROOT'};
        for iField=1:length(field_symbol)
            if ~isfield(x.D,field_symbol{iField})
                warning(['x.D.' field_symbol{iField} ' was not defined in DATA_PAR.m!'])
            end
        end
    end

    if ~isfield(x,'dataset')
        warning('x.dataset didn''nt exist');
    else
        if ~isempty(regexp(x.dataset.subjectRegexp, '^(\^|)\.\*(\$|)$', 'once'))
            warning('Subject regexp not specific! Check that no wrong folders are included as subjects');
        end
    end

    fprintf('\n');

end


