function EPAD_ASL_parmsPrepare(AnalysisDir)
%EPAD_ASL_parmsPrepare Prepare ASL parameter files for different EPAD vendors/sites
%
% FORMAT: EPAD_ASL_parmsPrepare(AnalysisDir)
% 
% INPUT:
%   AnalysisDir - path to folder containing the NIfTI/BIDS data (e.g. /data/RAD/share/EPAD/analysis) (REQUIRED)
%
% OUTPUT: n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function sets the specific EPAD ASL vendor parameters &
%              puts them in '_parms.mat' sidecars. REPLACE THIS FUNCTION BY EPAD_CREATEASLJSONPARS.M
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: EPAD_ASL_parmsPrepare('data/RAD/share/EPAD500/analysis);
% __________________________________
% Copyright 2015-2019 ExploreASL


SequenceName = {'Philips2DEPI1'     'Philips2DEPI2' 'Siemens2DEPI' 'Siemens3DGRASE'};
Sites        = {'(020|022|030|040)' '050'           '(012|031)'    '(010|011|060)'};

G.M0 				  = {'separate_scan' 'separate_scan' 'no_background_suppression' 'separate_scan'};
G.BackGrSupprPulses   = {2 2 0 2};
G.readout_dim         = {'2D' '2D' '2D' '3D'};
G.Vendor        	  = {'Philips' 'Philips' 'Siemens' 'Siemens'};
G.LabelingType        = {'CASL' 'CASL' 'PASL' 'PASL'};
G.Initial_PLD         = {2025 1800 2000 2000};
G.LabelingDuration    = {1650 1650 800 800};
G.SliceReadoutTime    = {36.5278 36.5278 35 0};

FieldsG = fields(G);

fprintf('%s','Preparing ASL parameter files for different sites:   ');
SubjectList = xASL_adm_GetFsList(AnalysisDir, '^\d{3}EPAD\d*(|_\d*)$', true, [], [], [0 Inf]);

if isempty(SubjectList)
    warning('No subjects found for ASL site-specific curation, skipping');
    return;
end

for iSubject=1:length(SubjectList)
    xASL_TrackProgress(iSubject,length(SubjectList));
    % Find which sequence the site uses
    SequenceN = find(cellfun(@(x) ~isempty(regexp(SubjectList{iSubject}(1:3),x)), Sites));
    if  length(SequenceN)==0
        error(['No sequence defined for ' SubjectList{iSubject}]);
    elseif length(SequenceN)>1
        error(['Multiple sequence defined for ' SubjectList{iSubject} ', should be only 1']);
    end
    
    % Find ASL&M0 JSON files
    JsonFiles = xASL_adm_GetFileList(fullfile(AnalysisDir, SubjectList{iSubject}), '^(ASL4D|M0)\.json$', 'FPListRec', [0 Inf]);
    for iJson=1:length(JsonFiles)    
        parms = xASL_import_json(JsonFiles{iJson});
        for iG=1:length(FieldsG)
            parms.(FieldsG{iG}) = G.(FieldsG{iG}){SequenceN};
        end
        xASL_adm_SaveJSON(parms, JsonFiles{iJson});
    end
    
end

xASL_TrackProgress(1, 1);
fprintf('\n');


end

