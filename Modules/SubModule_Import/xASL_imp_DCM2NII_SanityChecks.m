function xASL_imp_DCM2NII_SanityChecks(x)
%xASL_imp_DCM2NII_SanityChecks Sanity check for missing elements
%
% FORMAT: xASL_imp_DCM2NII_SanityChecks(x)
%
% INPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Sanity check for missing elements.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
%
% __________________________________
% Copyright 2015-2021 ExploreASL


    if x.modules.import.numOf.nSubjects==0
        error('No subjects')
    end
    if x.modules.import.numOf.nVisits==0
        error('No visits')
    end
    if x.modules.import.numOf.nSessions==0
        error('No sessions')
    end
    if x.modules.import.numOf.nScans==0
        error('No scans')
    end

end


