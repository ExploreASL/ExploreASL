function EPAD_CreateASLJSONPars(AnalysisDir)
%EPAD_ASL_parmsPrepare Prepare ASL parameter files for different EPAD vendors/sites
%
% FORMAT: EPAD_CreateASLJSONPars(AnalysisDir)
% 
% INPUT:
%   AnalysisDir - path to folder containing the NIfTI/BIDS data (e.g. /data/RAD/share/EPAD/analysis) (REQUIRED)
%
% INPUT file:
%   '//ExploreASL/CustomScripts/EPAD/DataPars/ASLQParms.json'
%   this file contains the "MotherDatabase" for the ASL vendor settings of
%   EPAD
%
% OUTPUT: n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function sets the specific EPAD ASL vendor parameters &
%              puts them in JSON sidecars
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: EPAD_CreateASLJSONPars('data/RAD/share/EPAD500/analysis);
% __________________________________
% Copyright 2015-2019 ExploreASL

if exist(fullfile(pwd,'ExploreASL_Master.m')) % we are in the ExploreASL root folder
    QParmsPath = fullfile('CustomScripts','EPAD','DataParsOther','ASLQParms.json');
elseif exist(fullfile(pwd,'EPAD_CreateASLJSONPars.m')) % we are in the EPAD CustomScripts folder
    QParmsPath = fullfile('DataParsOther','ASLQParms.json');
else
    warning('Couldnt find the mother database ASLQParms.json, skipping...');
    return;
end

MotherData = spm_jsonread(QParmsPath);
SubjectDirs = xASL_adm_GetFileList(AnalysisDir, '^\d{3}EPAD\d*$', 'FPList', [0 Inf], true);

fprintf('Populating ASL JSON files with vendor/sequence-specific quantification parameters:   ');

for iS=1:length(SubjectDirs)
    xASL_TrackProgress(iS,length(SubjectDirs));
    % Identify site
    [~, CurrentID] = fileparts(SubjectDirs{iS});
    CurrentSite = ['Site' CurrentID(1:3)];
    % Search for ASL JSON
    JSONPath = xASL_adm_GetFileList(fullfile(SubjectDirs{iS}, 'ASL_1'), '^ASL(?!.*RevPE).*\.json$', 'FPList', [0 Inf]);

    if ~isempty(JSONPath)
        for iC=1:length(JSONPath)
            jsonData = spm_jsonread(JSONPath{iC});
            FieldsAre = fields(MotherData.(CurrentSite));
            for iField=1:length(FieldsAre)
                jsonData.(FieldsAre{iField}) = MotherData.(CurrentSite).(FieldsAre{iField});
            end
            xASL_adm_SaveJSON(jsonData,JSONPath{iC});
        end
    end
end
    
fprintf('\n');
    
end



