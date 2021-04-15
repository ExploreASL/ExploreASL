function [dcm2niiCatchedErrors] = xASL_imp_CatchErrors(WarningID, WarningMessage, WarningLine, WarningFileName, WarningPath, scan_name, scanpath, destdir, dcm2niiCatchedErrors, imPar, StackIn)
%xASL_imp_CatchErrors Catch Errors.
%
% FORMAT: [dcm2niiCatchedErrors] = xASL_imp_CatchErrors(WarningID, WarningMessage, WarningLine, WarningFileName, WarningPath, scan_name, scanpath, destdir, dcm2niiCatchedErrors, imPar, StackIn)
% 
% INPUT:
%   WarningID            - Warning ID
%   WarningMessage       - Warning message
%   WarningLine          - Warning line
%   WarningFileName      - Warning file name
%   WarningPath          - Warning path
%   scan_name            - Scan name
%   scanpath             - scan path
%   destdir              - Destination directory
%   dcm2niiCatchedErrors - DCM2NII catched errors
%   imPar                - imPar struct
%   StackIn              - Stack In
%
% OUTPUT:
%   dcm2niiCatchedErrors - DCM2NII catched errors
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Catch reported warnings/errors, print them if verbose, & add them to a structure of warnings/errors to be stored for later QC.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     [dcm2niiCatchedErrors] = xASL_imp_CatchErrors(WarningID, WarningMessage, WarningLine, WarningFileName, WarningPath, scan_name, scanpath, destdir, dcm2niiCatchedErrors, imPar, StackIn);
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% Catch Errors

    if imPar.bVerbose % print warning if we want verbose
        warning(WarningMessage);
    end

    % Find index of the warning to store
    if isempty(fields(dcm2niiCatchedErrors))
        IndexN = 1;
    else
        IndexN = length(dcm2niiCatchedErrors)+1;
    end

    % store the warning/error
    dcm2niiCatchedErrors(IndexN).scan_name = scan_name;
    dcm2niiCatchedErrors(IndexN).scanpath = scanpath;
    dcm2niiCatchedErrors(IndexN).destdir = destdir;
    dcm2niiCatchedErrors(IndexN).identifier = WarningID;
    dcm2niiCatchedErrors(IndexN).message = WarningMessage;
    dcm2niiCatchedErrors(IndexN).cause = 'n/a';

    if exist('StackIn', 'var')
        dcm2niiCatchedErrors(IndexN).stack = StackIn;
    else
        dcm2niiCatchedErrors(IndexN).stack.file = fullfile(WarningPath, [WarningFileName '.m']);
        dcm2niiCatchedErrors(IndexN).stack.name = WarningFileName;
        dcm2niiCatchedErrors(IndexN).stack.line = WarningLine(end).line;
    end

end



