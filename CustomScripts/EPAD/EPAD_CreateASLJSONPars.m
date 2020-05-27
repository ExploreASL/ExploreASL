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



if exist(fullfile(pwd,'ExploreASL_Master.m'), 'file') % we are in the ExploreASL root folder
    QParmsPath = fullfile('CustomScripts','EPAD','DataParsOther','ASLQParms.json');
elseif exist(fullfile(pwd,'EPAD_CreateASLJSONPars.m'), 'file') % we are in the EPAD CustomScripts folder
    QParmsPath = fullfile('DataParsOther','ASLQParms.json');
else
    MyPath = which('ExploreASL_Master.m');
    if ~isempty(MyPath) && exist(MyPath, 'file')
        QParmsPath = fullfile(fileparts(MyPath), 'CustomScripts','EPAD','DataParsOther','ASLQParms.json');
    else
        QParmsPath = '';
    end
end

if ~exist(QParmsPath, 'file')
    warning('Couldnt find the mother database ASLQParms.json, skipping...');
    return;
end

MotherData = xASL_import_json(QParmsPath); % spm_jsonread
SubjectList = xASL_adm_GetFileList(AnalysisDir, '^\d{3}EPAD\d*(|_\d*)$', 'FPList', [0 Inf], true);

if isempty(SubjectList)
    warning('No subjects found for curating ASL metadata, skipping');
    return;
end

fprintf('Populating ASL JSON files with vendor/sequence-specific quantification parameters:   ');

for iSubject=1:length(SubjectList)
    xASL_TrackProgress(iSubject,length(SubjectList));
    % Identify site
    [~, CurrentID] = fileparts(SubjectList{iSubject});
    CurrentSite = ['Site' CurrentID(1:3)];
    % Search for ASL JSON
    JSONPath = xASL_adm_GetFileList(fullfile(SubjectList{iSubject}, 'ASL_1'), '^ASL(?!.*RevPE).*\.json$', 'FPList', [0 Inf]);

    if ~isempty(JSONPath)
        for iJSON=1:length(JSONPath)
            jsonData = spm_jsonread(JSONPath{iJSON});
            FieldsAre = fields(MotherData.(CurrentSite));
            for iField=1:length(FieldsAre)
                jsonData.(FieldsAre{iField}) = MotherData.(CurrentSite).(FieldsAre{iField});
            end
            xASL_adm_SaveJSON(jsonData,JSONPath{iJSON});
        end
    end
end
   
xASL_TrackProgress(1, 1);
fprintf('\n');
    
end



