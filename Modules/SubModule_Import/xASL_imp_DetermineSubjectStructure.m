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



    if x.opts.ImportModules(1)
        %% Import sourcedata structure
        if isfield(x.modules.import,'imPar') && isstruct(x.modules.import.imPar)
            % Basic import checks before execution
            x = xASL_imp_CheckImportSettings(x);
            
            % Check directories and permissions
            x = xASL_imp_CheckDirectoriesAndPermissions(x);
            
            % Here we try to fix backwards compatibility
            x.modules.import.imPar = xASL_imp_TokenBackwardsCompatibility(x.modules.import.imPar);

            % Read sourcedata
            x = xASL_imp_ReadSourceData(x);
            
            % Determine structure from sourcedata
            x = xASL_imp_DetermineStructureFromSourcedata(x);
            
            % Sanity check for missing elements
            xASL_imp_DCM2NII_SanityChecks(x);
            
            % Preallocate space for (global) counts
            x = xASL_imp_PreallocateGlobalCounts(x);
        else
            error('The imPar struct does not exist...');
        end
    elseif x.opts.ImportModules(2)
        %% Import temp data structure
        if isfield(x.modules.import,'imPar') && isstruct(x.modules.import.imPar)
            % Basic import checks before execution
            x = xASL_imp_CheckImportSettings(x);
            
            % Check directories and permissions
            x = xASL_imp_CheckDirectoriesAndPermissions(x);
            
            % Determine structure from temp data
            x = xASL_imp_DetermineStructureFromTempdata(x);
        else
            error('The imPar struct does not exist...');
        end
    elseif x.opts.ImportModules(3) || x.opts.ImportModules(4)
        %% Import the rawdata structure
        if isfield(x.modules.import,'imPar') && isstruct(x.modules.import.imPar)
            % Basic import checks before execution
            x = xASL_imp_CheckImportSettings(x);
            
            % Check directories and permissions
            x = xASL_imp_CheckDirectoriesAndPermissions(x);
            
            % Determine structure from temp data
            x = xASL_imp_DetermineStructureFromRawdata(x);
        else
            error('The imPar struct does not exist...');
        end
    else
        %% Not a valid option
        error('Invalid option of import settings...');
    end



end


