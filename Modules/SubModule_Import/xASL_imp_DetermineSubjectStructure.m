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
        % Import sourcedata structure
        if ~isnan(x.modules.import.imPar)

        else
            error('The imPar struct does not exist...');
        end
    elseif x.opts.ImportModules(2)
        % Import temp data structure
        if ~isnan(x.modules.import.imPar)

        else
            error('The imPar struct does not exist...');
        end
    else
        % Not a valid option
        error('Invalid option of import settings...');
    end



end



