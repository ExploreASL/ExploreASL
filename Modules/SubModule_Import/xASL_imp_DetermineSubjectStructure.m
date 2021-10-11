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
        % Determine structure from temp data
        x = xASL_imp_DetermineStructureFromRawdata(x);
        
    end
    
    % SESSIONS DUMMY
    % fprintf(2,'Currently x.SESSIONS is not supported for the import...\n');
    x.SESSIONS = {''};


end




