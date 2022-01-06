function xASL_imp_DCM2NII_SanityChecks(x,thisSubject,thisVisit)
%xASL_imp_DCM2NII_SanityChecks Sanity check for missing elements
%
% FORMAT: xASL_imp_DCM2NII_SanityChecks(x)
%
% INPUT:
%   x           - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%   thisSubject - Current subject x.overview.(sFieldName)
%   thisVisit   - Current visit x.overview.(sFieldName).(vFieldName)
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


    if x.modules.import.nSubjects==0
        error('No subjects')
    end
    if thisSubject.nVisits==0
        error('No visits')
    end
    if thisVisit.nSessions==0
        error('No sessions')
    end
    if thisVisit.nScans==0
        error('No scans')
    end

end


