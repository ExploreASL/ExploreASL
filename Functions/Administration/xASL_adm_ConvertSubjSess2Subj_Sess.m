function [iSubj iSess] = xASL_adm_ConvertSubjSess2Subj_Sess(nSessions, iSubjSess)
%xASL_adm_ConvertSubjSess2Subj_Sess Converts combined SubjectSession index to subject &
%session indices. Useful for data lists in ExploreASL
%
% FORMAT:       [iSubj iSess] = xASL_adm_ConvertSubjSess2Subj_Sess(nSessions, iSubjSess)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Converts combined SubjectSession index to subject & session
%               indices. Useful for data lists in ExploreASL.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

    iSubj    = ceil(iSubjSess/nSessions);
    iSess    = iSubjSess - (iSubj-1)*nSessions;
end
