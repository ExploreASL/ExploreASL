function x = xASL_imp_PreallocateGlobalCounts(x)
%xASL_imp_PreallocateGlobalCounts Preallocate space for (global) counts
%
% FORMAT: x = xASL_imp_PreallocateGlobalCounts(x)
%
% INPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Preallocate space for (global) counts.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
%
% __________________________________
% Copyright 2015-2021 ExploreASL


    % keep a count of all individual scans
    x.modules.import.globalCounts.converted_scans = ...
        zeros(x.modules.import.numOf.nSubjects, x.modules.import.numOf.nVisits, x.modules.import.numOf.nSessions, x.modules.import.numOf.nScans,'uint8');
    % keep a count of all individual scans
    x.modules.import.globalCounts.skipped_scans = ...
        zeros(x.modules.import.numOf.nSubjects, x.modules.import.numOf.nVisits, x.modules.import.numOf.nSessions, x.modules.import.numOf.nScans,'uint8');
    % keep a count of all individual scans
    x.modules.import.globalCounts.missing_scans = ...
        zeros(x.modules.import.numOf.nSubjects, x.modules.import.numOf.nVisits, x.modules.import.numOf.nSessions, x.modules.import.numOf.nScans,'uint8');
    
    % define a cell array for storing info for parameter summary file
    x.modules.import.summary_lines = ...
        cell(x.modules.import.numOf.nSubjects, x.modules.import.numOf.nVisits, x.modules.import.numOf.nSessions, x.modules.import.numOf.nScans);


end



