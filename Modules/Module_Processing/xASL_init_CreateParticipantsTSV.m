function [x] = xASL_init_CreateParticipantsTSV(x)
%xASL_init_CreateParticipantsTSV Create participants.tsv file
%
% FORMAT: [x] = xASL_init_CreateParticipantsTSV(x)
%
% INPUT:
%   x       - ExploreASL x structure (STRUCT, REQUIRED)
%
% OUTPUT:
%   x       - ExploreASL x structure (STRUCT)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Load participants.tsv file (which contains metadata for the BIDS dataset
%
% EXAMPLE:     This is part of the initialization workflow. Check out the usage there.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:  n/a
% __________________________________
% Copyright (c) 2015-2024 ExploreASL



    %% Create participants.tsv, if needed
    
    % Define file names
    pathParticipantsTSV = fullfile(x.dir.DatasetRoot, 'participants.tsv');
    pathParticipantsTSVxASL = fullfile(x.dir.xASLDerivatives, 'participants.tsv');
    
    % Backwards compatibility: rename to lowercase
    fileListParticipantsTSVold = xASL_adm_GetFileList(x.dir.DatasetRoot,'^Participants.tsv$',false);
    if ~isempty(fileListParticipantsTSVold) % We added the _temp copy step so that the code works on case insensitive systems like windows as well. Please don't remove that step for backwards compatibility (at least not until release 2.0.0).
        pathParticipantsTSVold = fullfile(x.dir.DatasetRoot, 'Participants.tsv');
        pathParticipantsTSVoldTemp = fullfile(x.dir.DatasetRoot, 'Participants_temp.tsv');
        xASL_Move(pathParticipantsTSVold,pathParticipantsTSVoldTemp);
        xASL_Move(pathParticipantsTSVoldTemp,pathParticipantsTSV);
    end
    
    % Check if participants.tsv exists & copy it to the derivatives
    if ~xASL_exist(pathParticipantsTSVxASL, 'file')
        xASL_Copy(pathParticipantsTSV, pathParticipantsTSVxASL);
    else
        warning([pathParticipantsTSVxASL ' already existed, not overwritten']);
    end


end