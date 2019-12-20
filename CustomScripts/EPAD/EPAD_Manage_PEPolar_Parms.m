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
%              ExploreASL_Import). We assume that these are better than
%              those inside the JSON, hence replace them. Here the phase
%              encoding parameters are added specified per vendor.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: EPAD_Manage_PEPolar_Parms(AnalysisDir);
% __________________________________
% Copyright 2015-2019 ExploreASL


PEPolarParms = {'EffectiveEchoSpacing' 'AcquisitionMatrix' 'TotalReadoutTime'};

SubjectList = xASL_adm_GetFsList(AnalysisDir,'^\d{3}EPAD\d*$', true, [], [], [0 Inf]);

fprintf('%s','Managing PEPolar parameters:  0%');

for iSubject=1:length(SubjectList)
    xASL_TrackProgress(iSubject,length(SubjectList));
    CurrDir = fullfile(AnalysisDir, SubjectList{iSubject});
    Filelist = xASL_adm_GetFileList(CurrDir, '.*_parms\.mat$', 'FPListRec', [0 Inf], false);
    
    for iFile=1:length(Filelist)
        %% First check if a JSON sidecar exists:
        ParmsPath = Filelist{iFile};
        IndexN = length('_parms.mat');
        JSONPath = [ParmsPath(1:end-IndexN) '.json'];
        
        Ind1 = regexp(JSONPath,'run(-|_)\d*');
        if ~exist(JSONPath,'file') && ~isempty(Ind1)
            JSONPath(Ind1+3) = '_';
            if ~exist(JSONPath,'file')
                JSONPath(Ind1+3) = '-';
                if ~exist(JSONPath,'file')
                    continue;
                end
            end
        end
        
        if exist(JSONPath,'file') && exist(ParmsPath,'file')
            %% Load parms.mat
            Parms = load(ParmsPath,'-mat');
            
            % clear any previous fields
            IsComplete = true;
            for iField=1:length(PEPolarParms)
                if isfield(Parms,PEPolarParms{iField})
                    % First remove the field if it already exists
                    Parms = rmfield(Parms,PEPolarParms{iField});
                end
                if isfield(Parms.parms,PEPolarParms{iField})
                    % Now copy the field
                    Parms.(PEPolarParms{iField}) = Parms.parms.(PEPolarParms{iField});
                    if ~isfinite(Parms.(PEPolarParms{iField}))
                        IsComplete = false;
                    end
                else
                    Parms.(PEPolarParms{iField}) = NaN;
                    IsComplete = false;
                end
            end
            
            if IsComplete
                %% Load the JSON
                jsonData = spm_jsonread(JSONPath); % read the JSON
                
                for iField=1:length(PEPolarParms)
                    jsonData.(PEPolarParms{iField}) = Parms.(PEPolarParms{iField});
                end
               
                if ~isfield(jsonData,'PhaseEncodingDirection')
                    if strcmp(jsonData.Manufacturer,'Philips')
                        if ~isempty(strfind(lower(JSONPath),'asl')) || ~isempty(strfind(lower(JSONPath),'func'))
                            if ~isempty(strfind(lower(JSONPath),'revpe'))
                                jsonData.PhaseEncodingDirection = 'j-'; % for EPAD ASL & func RevPE
                            else
                                jsonData.PhaseEncodingDirection = 'j'; % for EPAD ASL & func NormPE
                            end
                        elseif ~isempty(strfind(lower(JSONPath),'dwi'))
                            if ~isempty(strfind(lower(JSONPath),'revpe'))
                                jsonData.PhaseEncodingDirection = 'j'; % for EPAD DTI RevPE
                            else
                                jsonData.PhaseEncodingDirection = 'j-'; % for EPAD DTI NormPE
                            end
                        end
                    end
                end
                % Siemens should have the same, but these PhaseEncodingDirection were populated by dcm2niiX     

                % NB: NEED TO ADD THE FMRI DIRECTION

                fclose all;
                xASL_delete(JSONPath);
                xASL_adm_SaveJSON(jsonData,JSONPath); % write/update the JSON
            end
        end
    end
end

fprintf('\n');

end