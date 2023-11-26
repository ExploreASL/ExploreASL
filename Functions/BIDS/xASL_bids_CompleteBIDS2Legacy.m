function x = xASL_bids_CompleteBIDS2Legacy(x)
%xASL_bids_CompleteBIDS2Legacy Finish the remaining conversion steps, which are not done for each subject individually
%
% FORMAT: x = xASL_bids_CompleteBIDS2Legacy(x)
%
% INPUT:
%   x      - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x      - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Finish the remaining conversion steps, which are not done for each subject individually.
%
% 2. Copy participants.tsv
% 3. Copy dataset_description JSON
% 4. Add missing fields
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        x = xASL_bids_CompleteBIDS2Legacy(x);
%
% __________________________________
% Copyright 2015-2023 ExploreASL



%% Create participants.tsv

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
if xASL_exist(pathParticipantsTSV,'file')
    xASL_Copy(pathParticipantsTSV,pathParticipantsTSVxASL);
end



%% Create dataset_description.json

% Copy dataset_description JSON file
xASL_Copy(fullfile(x.dir.RawData, 'dataset_description.json'), fullfile(x.dir.xASLDerivatives, 'dataset_description.json'));


end