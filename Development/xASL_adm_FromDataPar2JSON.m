function xASL_adm_FromDataPar2JSON(DataParPath)
%xASL_adm_FromDataPar2JSON This function takes all parameters from the
%DataPar & moves them into all lower-level JSONs, per BIDS inheritance
%Note that this function assumes that the DataPar file is in the ROOT
%folder of the study, that contains all the JSON sidecars

%% 1) Load DataPar file
if nargin<1 || isempty(DataParPath) || ~exist(DataParPath, 'file')
    error('Invalid input argument or non-existing DataPar file');
end

DataPar = xASL_import_json(DataParPath);
FieldsAre = fields(DataPar);

%% 2) Get list of JSONs
AnalysisDir = fileparts(DataParPath);

fprintf('Obtaining list of JSON files\n');
FileList = xASL_adm_GetFileList(AnalysisDir, '^.*\.json$','FPListRec',[0 Inf]);

fprintf('Processing JSON files:   \n');

%% 3) Process the JSONs
for iFile=1:length(FileList)
    xASL_TrackProgress(iFile, length(FileList));
    [Fpath, Ffile] = xASL_fileparts(FileList{iFile});
    niiPath = fullfile(Fpath, [Ffile '.nii']);
    if ~xASL_exist(niiPath, 'file')
        % skip this JSON, it is not a sidecar
    else
        JSON = xASL_import_json(FileList{iFile});
        JSON = InsertFields(DataPar, FieldsAre, JSON);
        
        spm_jsonwrite(FileList{iFile}, JSON);
    end
end

fprintf('\n');

end


function JSON = InsertFields(DataPar, FieldsAre, JSON)
%InsertFields

Fields2Skip = {'Quality' 'DELETETEMP' 'subject_regexp' 'name'};
% These fields are environment parameters, not ASL-specific parameters

    for iField=1:length(FieldsAre)
        FieldValue = DataPar.(FieldsAre{iField});
        if ischar(FieldValue) || isnumeric(FieldValue) || islogical(FieldValue)
            SkipField = max(cellfun(@(y) strcmp(FieldsAre{iField},y), Fields2Skip));
            
            if isfield(JSON,FieldsAre{iField}) || SkipField
                % Skip this field: per inheritance principle, daughters
                % fields have preference
                % Also skip this field if it is an environment parameter
            else
                JSON.(FieldsAre{iField}) = DataPar.(FieldsAre{iField});
            end
        else % assume we have subfields
            try
                Subfields = fields(DataPar.(FieldsAre{iField}));
            catch ME
                warning(['Something went wrong with field ' FieldsAre{iField}]);
                fprintf('Is this field not char/numeric and it doesnt have subfields?\n');
                fprintf('%s\n', ME.message);
                continue; % with next field
            end
            % do the same for the subfields
            if ~isfield(JSON,FieldsAre{iField})
                JSON.(FieldsAre{iField}) = struct;
            end
            
            JSON.(FieldsAre{iField}) = InsertFields(DataPar.(FieldsAre{iField}), Subfields, JSON.(FieldsAre{iField}));
        end
    end


end

