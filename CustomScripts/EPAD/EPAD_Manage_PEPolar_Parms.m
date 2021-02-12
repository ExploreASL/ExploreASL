function EPAD_Manage_PEPolar_Parms(AnalysisDir)
%EPAD_Manage_PEPolar_Parms Manage the PE parameters for TopUp
%
% FORMAT: EPAD_Manage_PEPolar_Parms(AnalysisDir)
% 
% INPUT:
%   AnalysisDir - path to folder containing the NIfTI/BIDS data (e.g. /data/RAD/share/EPAD/analysis) (REQUIRED)
%
% OUTPUT: n/a
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function takes the phase encoding (PE) parameters for
%              TopUp from the '_parms.mat' file (i.e. the suffix of the
%              sidecar file that should accompany the NIfTI after
%              ExploreASL_Import) if there is one. We assume that these are better than
%              those inside the JSON, hence replace them. 
%              Also, here the vndor-specific phase
%              encoding parameters are added.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: EPAD_Manage_PEPolar_Parms(AnalysisDir);
% __________________________________
% Copyright 2015-2019 ExploreASL


PEPolarParms = {'EffectiveEchoSpacing' 'AcquisitionMatrix' 'TotalReadoutTime'};

SubjectList = xASL_adm_GetFsList(AnalysisDir,'^\d{3}EPAD\d*(|_\d*)$', true, [], [], [0 Inf]);

if isempty(SubjectList)
    warning('Didnt find subjects for PEPolar curation, skipping');
    return;
end

fprintf('%s','Managing PEPolar parameters:  0%');

for iSubject=1:length(SubjectList)
    xASL_TrackProgress(iSubject,length(SubjectList));
    CurrDir = fullfile(AnalysisDir, SubjectList{iSubject});
    Filelist = xASL_adm_GetFileList(CurrDir, '^(?!(swi|T1|FLAIR|anat)).*\.json$', 'FPListRec', [0 Inf], false); % here we should search for JSON instead
    
    for iFile=1:length(Filelist)
        SaveJSON = false;
        %% First check if a parms.mat sidecar exists, then let this overwrite the JSON
        JSONPath = Filelist{iFile};
        [Fpath, jsonFile] = fileparts(JSONPath);
        ParmsPath = fullfile(Fpath, [jsonFile '_parms.mat']);
        jsonData = xASL_import_json(JSONPath); % read the JSON
        
        %% Load parms.mat
        if exist(ParmsPath,'file')
            Parms = load(ParmsPath, '-mat');
            for iField=1:length(PEPolarParms)
                % if we have meaningfull information in the _parms.mat sidecar
                if isfield(Parms, PEPolarParms{iField}) && isfinite(Parms.PEPolarParms{iField})
                    % overwrite the JSON
                    jsonData.(PEPolarParms{iField}) = Parms.parms.(PEPolarParms{iField});
                    SaveJSON = true;
                end
            end
        end
        
        %% Deal with Philips defaults in EPAD
        if ~isfield(jsonData, 'Manufacturer')
            warning('Couldnt find Manufacturer field, skipping curating TopUp parameters');
            continue;
        elseif ~strcmp(jsonData.Manufacturer, 'Philips')
            if ~isfield(jsonData, 'PhaseEncodingDirection')
                warning('PhaseEncodingDirection field missing for:');
                fprintf('%s\n', ParmsPath);
            end            
            
            % Skip this, below is only for Philips
            % Siemens information is correctly found in the DICOMs by
            % dcm2niiX, GE we dont have in EPAD
            continue;
        end

        if ~isempty(strfind(lower(jsonFile), 'asl')) || ~isempty(strfind(lower(jsonFile), 'm0')) || ~isempty(strfind(lower(jsonFile), 'func')
            if ~isempty(strfind(lower(jsonFile),'revpe'))
                jsonData.PhaseEncodingDirection = 'j-'; % for EPAD ASL & func RevPE
            else
                jsonData.PhaseEncodingDirection = 'j'; % for EPAD ASL & func NormPE
            end
            SaveJSON = true;
        elseif ~isempty(strfind(lower(jsonFile),'dwi'))
            if ~isempty(strfind(lower(jsonFile),'revpe'))
                jsonData.PhaseEncodingDirection = 'j'; % for EPAD DTI RevPE
            else
                jsonData.PhaseEncodingDirection = 'j-'; % for EPAD DTI NormPE
            end
            SaveJSON = true;
        end

        if SaveJSON
            xASL_delete(JSONPath);
            xASL_adm_SaveJSON(jsonData, JSONPath); % write/update the JSON
        end
    end
end

xASL_TrackProgress(1, 1);
fprintf('\n');

end