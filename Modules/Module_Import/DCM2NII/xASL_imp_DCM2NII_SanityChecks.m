function xASL_imp_DCM2NII_SanityChecks(x, thisSubject, thisVisit)
%xASL_imp_DCM2NII_SanityChecks Sanity check for missing elements
%
% FORMAT: xASL_imp_DCM2NII_SanityChecks(x, thisSubject, thisVisit)
%
% INPUT:
%   x           - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%   thisSubject - Current subject x.importOverview.(sFieldName)
%   thisVisit   - Current visit x.importOverview.(sFieldName).(vFieldName)
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
% Copyright 2015-2022 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________



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


