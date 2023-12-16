function [subject,visit] = xASL_imp_PreallocateGlobalCounts(nSubjects, subject, visit)
%xASL_imp_PreallocateGlobalCounts Preallocate space for (global) counts
%
% FORMAT: x = xASL_imp_PreallocateGlobalCounts(nSubjects, subject, visit)
%
% INPUT:
%   nSubjects - Number of subjects (INTEGER)
%   subject   - Current subject x.importOverview.(sFieldName)
%   visit     - Current visit x.importOverview.(sFieldName).(vFieldName)
%
% OUTPUT:
%   visit - Current visit x.importOverview.(sFieldName).(vFieldName)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Preallocate space for (global) counts.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
%
% __________________________________
% Copyright 2015-2022 ExploreASL


    % keep a count of all individual scans
    subject.globalCounts.converted_scans = zeros(nSubjects, subject.nVisits, visit.nSessions, visit.nScans,'uint8');
    
    % keep a count of all individual scans
    subject.globalCounts.skipped_scans = zeros(nSubjects, subject.nVisits, visit.nSessions, visit.nScans,'uint8');
    
    % keep a count of all individual scans
    subject.globalCounts.missing_scans = zeros(nSubjects, subject.nVisits, visit.nSessions, visit.nScans,'uint8');
    
    % define a cell array for storing info for parameter summary file
    subject.summary_lines = cell(nSubjects, subject.nVisits, visit.nSessions, visit.nScans);


end



