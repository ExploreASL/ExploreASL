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



    iSubj    = ceil(iSubjSess/nSessions);
    iSess    = iSubjSess - (iSubj-1)*nSessions;

end

