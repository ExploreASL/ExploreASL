function xASL_bids_GenerateParticipantsJSON(x)
%xASL_bids_GenerateParticipantsJSON This function generates the participants.json sidecar to the participants.tsv
%
% FORMAT: xASL_bids_GenerateParticipantsJSON(x)
%
% INPUT:
%   x       - x structure containing all input parameters (REQUIRED)
%   participants.tsv needs to be present
%
% OUTPUT: n/a 
%   participants.json will be written
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function generates the participants.json sidecar to the participants.tsv, using the following steps
%
% 1. Load pre-existing participants.json
% 2. Parse participants.tsv
% 3. Parse standard keys
% 4. Parse ASL-related keys
% 4a. Replace key abbreviations with verbose description
% 4b. Manage units
% 5. (over-)Write the new json
%
% EXAMPLE: xASL_bids_GenerateParticipantsJSON(x);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2025 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________



%% 1. Load pre-existing participants.json
pathJson = fullfile(x.dir.xASLDerivatives, 'participants.json');
pathTSV = fullfile(x.dir.xASLDerivatives, 'participants.tsv');
if exist(pathJson, 'file')
    % Load the JSON, but note that we will overwrite it if fields already exist
    jsonParticipants = xASL_io_ReadJson(pathJson);
else
    jsonParticipants = struct();
end


%% 2. Parse participants.tsv
if ~exist(pathTSV, 'file')
    warning(['Skipping, file missing: ' pathTSV]);
    return;
end

tableParticipants = xASL_tsvRead(pathTSV);
listKeys = tableParticipants(1,:);

for iList=1:length(listKeys)
    switch lower(listKeys{iList})

        
        %% 3. Parse standard keys
        case 'participant_id'
            jsonParticipants.participant_id.Description = 'Unique participant identifier, including _session suffix (==visit in ExploreASL legacy)';
        case 'session'
            jsonParticipants.session.Description = 'Unique run identifier (== session in ExploreASL legacy)';
        case 'age'
            jsonParticipants.age.Description = 'Age of the participant at the time of scanning';
            jsonParticipants.age.Units = 'years';
        case 'sex'
            jsonParticipants.sex.Description = 'Biological sex of the participant';
            jsonParticipants.sex.Levels.F = 'female';
            jsonParticipants.sex.Levels.M = 'male';
        case 'site'
            jsonParticipants.site.Description = 'Center where this scan was made';
            uniqueSites = unique(tableParticipants(2:end, iList));
            for iSite=1:length(uniqueSites)
                jsonParticipants.site.Levels.(['Site_' num2str(iSite)]) = uniqueSites{iSite};
            end
        case 'gm_vol'
            jsonParticipants.gm_vol.Description = 'Gray matter volume of the participant at this visit';
            jsonParticipants.gm_vol.Units = 'liter';
        case 'wm_vol'
            jsonParticipants.wm_vol.Description = 'White matter volume of the participant at this visit';
            jsonParticipants.wm_vol.Units = 'liter';                    
        case 'csf_vol'
            jsonParticipants.csf_vol.Description = 'Cerebrospinal fluid volume of the participant at this visit';
            jsonParticipants.csf_vol.Units = 'liter';
        case {'gm_icv_ratio', 'gm_icvratio'}
            jsonParticipants.gm_icv_ratio.Description = 'Ratio of gray matter volume to the intracranial volume of the participant at this visit, proxy of gray matter atrophy';
            jsonParticipants.gm_icv_ratio.Units = 'ratio';
        case {'gmwm_icv_ratio', 'gmwm_icvratio'}
            jsonParticipants.gmwm_icv_ratio.Description = 'Ratio of parenchymal volume (gray matter + WM volume) to the intracranial volume of the participant at this visit, proxy of wholebrain atrophy';
            jsonParticipants.gmwm_icv_ratio.Units = 'ratio';
        case 'wmh_vol'
            jsonParticipants.wmh_vol.Description = 'White matter hyperintensity volume of the participant at this visit';
            jsonParticipants.wmh_vol.Units = 'liter';
        case 'wmh_count'
            jsonParticipants.wmh_count.Description = 'Number of non-connected white matter hyperintensities of the participant at this visit';
            jsonParticipants.wmh_count.Units = 'integer';
        case 'meanmotion'
            jsonParticipants.motion.Description = 'Mean net displacement vector difference (head motion) of this ASL scan';
            jsonParticipants.motion.Units = 'mm RMS';
        otherwise


            %% 4. Parse ASL-related keys
            if ~isempty(regexpi(listKeys{iList}, '(cbf|att|m0|control|tex|itt)'))
                descriptionBIDS = lower(listKeys{iList});


                %% 4a. Replace key abbreviations with verbose description
                descriptionBIDS = strrep(descriptionBIDS, 'qcbf', 'cerebral blood flow');
                descriptionBIDS = strrep(descriptionBIDS, 'cbf', 'cerebral blood flow');
                descriptionBIDS = strrep(descriptionBIDS, 'att', 'arterial transit time');
                descriptionBIDS = strrep(descriptionBIDS, 'm0', 'M0 reference intensity');
                descriptionBIDS = strrep(descriptionBIDS, 'meancontrol', 'control');
                descriptionBIDS = strrep(descriptionBIDS, 'itt', 'intravoxel transit time');
                descriptionBIDS = strrep(descriptionBIDS, 'tex', 'time of exchange');
                
                descriptionBIDS = strrep(descriptionBIDS, '_b', ' bilateral');
                descriptionBIDS = strrep(descriptionBIDS, '_l', ' left');
                descriptionBIDS = strrep(descriptionBIDS, '_r', ' right');
                descriptionBIDS = strrep(descriptionBIDS, '_aca', ' anterior cerebral artery territory');
                descriptionBIDS = strrep(descriptionBIDS, '_mca', ' middle cerebral artery territory');
                descriptionBIDS = strrep(descriptionBIDS, '_pca', ' posterior cerebral artery territory');
                descriptionBIDS = strrep(descriptionBIDS, '_gm', ' gray matter');
                descriptionBIDS = strrep(descriptionBIDS, '_wm', ' white matter');
                descriptionBIDS = strrep(descriptionBIDS, 'deepwm white matter', 'deep white matter');
                descriptionBIDS = strrep(descriptionBIDS, '_pvc0', ', without partial volume correction');
                descriptionBIDS = strrep(descriptionBIDS, '_pvc2', ', with partial volume correction');
                descriptionBIDS = strrep(descriptionBIDS, 'mean_', 'Mean ');
                descriptionBIDS = strrep(descriptionBIDS, 'cov_', 'Spatial coefficient of variation ');
                descriptionBIDS = strrep(descriptionBIDS, '_', ' ');

                

                %% 4b. Manage units
                if ~isempty(regexpi(descriptionBIDS, '(cbf|cerebral blood flow)'))
                    jsonParticipants.(lower(listKeys{iList})).Units = 'mL/100g/min';
                elseif ~isempty(regexpi(descriptionBIDS, ('abv|arterial blood volume)')))
                    jsonParticipants.(lower(listKeys{iList})).Units = '%';
                elseif ~isempty(regexpi(descriptionBIDS, '(att|tt|itt|tex|arterial transit time|intravoxel transit time|time of exchange)'))
                    jsonParticipants.(lower(listKeys{iList})).Units = 's';
                elseif ~isempty(regexpi(descriptionBIDS, '(m0|control)'))
                    jsonParticipants.(lower(listKeys{iList})).Units = 'a.u.';
                else
                    warning(['Unknown type, cannot establish the unit for ' listKeys{iList}]);
                    jsonParticipants.(lower(listKeys{iList})).Units = 'NaN';
                end

                jsonParticipants.(lower(listKeys{iList})).Description = descriptionBIDS;
            else
                warning(['Unknown participants.tsv field: ' listKeys{iList}]);
            end
    end
end


%% 5. (over-)Write the new json
xASL_io_WriteJson(pathJson, jsonParticipants);


end