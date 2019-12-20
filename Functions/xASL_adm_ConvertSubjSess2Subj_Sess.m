function [iSubj iSess] = xASL_adm_ConvertSubjSess2Subj_Sess(nSessions, iSubjSess)
%xASL_adm_ConvertSubjSess2Subj_Sess Converts combined SubjectSession index to subject &
%session indices. Useful for data lists in ExploreASL

    iSubj    = ceil(iSubjSess/nSessions);
    iSess    = iSubjSess - (iSubj-1)*nSessions;

end

