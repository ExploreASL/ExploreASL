function x = xASL_imp_FinishImport(x)
%xASL_imp_FinishImport Finish the remaining import steps, which are not done for each subject individually
%
% FORMAT: x = xASL_imp_FinishImport(x)
%
% INPUT:
%   x      - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x      - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Finish the remaining import steps, which are not done for each subject individually.
%
% 1. Create dataPar.json
% 2. Copy participants.tsv
% 3. Copy dataset_description JSON and add 'GeneratedBy' fields
% 4. Add missing fields
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        x = xASL_imp_FinishImport(x);
%
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% 1. Create dataPar.json
    if x.opts.ImportModules(4)
       [x,dataPar] = xASL_imp_FinishImport_CreateDataPar(x);
    end


    %% 2. Copy participants.tsv
    if x.opts.ImportModules(4)
        x = xASL_imp_FinishImport_CreateParticipants(x);
    end
    
    
    %% 3. Copy dataset_description JSON and add 'GeneratedBy' fields
    if x.opts.ImportModules(4)
        x = xASL_imp_FinishImport_CreateDatasetDescription(x);
    end
    
    
    %% 4. Add missing fields
    if x.opts.ImportModules(4)
        x = xASL_imp_FinishImport_AddMissingFields(x,dataPar);        
    end
    
    
end



%% Create the dataPar.json 
function [x,dataPar] = xASL_imp_FinishImport_CreateDataPar(x)

    if ~isfield(x,'dataPar')
        % Create default if no dataPar was provided
        fprintf('No dataPar.json provided, will create default version...\n');
        dataPar = struct();
    else
        % dataPar was provided
        fprintf('Export provided dataPar.json...\n');
        dataPar = x.dataPar;
        x = rmfield(x,'dataPar');
    end

    % Fills in important information in the dataPar if missing
    if ~isfield(dataPar,'x')
        % Add x field
        dataPar.x = struct;
    end
    % Dataset fields
    if ~isfield(dataPar.x,'dataset')
        dataPar.x.dataset = struct;
    end
    % Add default subject regular expression
    if ~isfield(dataPar.x.dataset,'subjectRegexp')
        dataPar.x.dataset.subjectRegexp = '^sub-.*$';
    end
    % Check for settings fields
    if ~isfield(dataPar.x,'settings')
        dataPar.x.settings = struct;
    end
    % Check for quality field
    if ~isfield(dataPar.x.settings,'Quality')
        dataPar.x.settings.Quality = 1;
    end
    % Check for DELETETEMP field
    if ~isfield(dataPar.x.settings,'DELETETEMP')
        dataPar.x.settings.DELETETEMP = 1;
    end

    % Write DataParFile if it does not exist already
    fListDataPar = xASL_adm_GetFileList(fullfile(x.dir.DatasetRoot,'derivatives','ExploreASL'),'(?i)(^dataPar.*\.json$)', 'FPList', [], 0);
    if isempty(fListDataPar)
        fprintf('Creating dataPar.json since file does not exist in derivatives directory...\n');
        spm_jsonwrite(fullfile(fullfile(x.dir.DatasetRoot,'derivatives','ExploreASL'), 'dataPar.json'), dataPar);
    end

    % Update dataPar path
    fListDataParLegacy = xASL_adm_GetFileList(fullfile(x.dir.DatasetRoot,'derivatives','ExploreASL'),'(?i)(^dataPar.*\.json$)', 'FPList', [], 0);
    x.dir.dataPar = fListDataParLegacy{1};
    if length(fListDataParLegacy)>1
        fprintf('Multiple dataPar.json files within the derivatives directory...\n');
    end

    % Overwrite dataPar.json in x structure
    fprintf('Overwriting x.dir.dataPar...\n');

end


%% Create participants.tsv
function x = xASL_imp_FinishImport_CreateParticipants(x)

    % Define file names
    pathParticipantsTSV = fullfile(x.dir.DatasetRoot, 'participants.tsv');
    pathParticipantsTSVxASL = fullfile(fullfile(x.dir.DatasetRoot,'derivatives','ExploreASL'), 'participants.tsv');

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

end


%% Create dataset_description.json
function x = xASL_imp_FinishImport_CreateDatasetDescription(x)

    try
        % Copy dataset_description JSON file
        xASL_Copy(fullfile(x.dir.DatasetRoot, 'rawdata', 'dataset_description.json'), ...
            fullfile(fullfile(x.dir.DatasetRoot,'derivatives','ExploreASL'), 'dataset_description.json'));

        % Get all JSON files in derivatives and add the "GeneratedBy" field
        jsonFiles = xASL_adm_GetFileList(fullfile(x.dir.DatasetRoot,'derivatives','ExploreASL'), '.json$', 'FPListRec');

        % Check file list
        if size(jsonFiles,1)>0
            % Iterate over files
            for iFile = 1:size(jsonFiles,1)
                % Check if file should be excluded first (exclude participants.json and dataPar.json)
                [~,fileName] = xASL_fileparts(jsonFiles{iFile});
                if isempty(regexpi(fileName,'participants')) && isempty(regexpi(fileName,'dataPar'))
                    thisPath = jsonFiles{iFile};
                    xASL_bids_AddGeneratedByField(x, thisPath);
                end
            end
        end
    catch ME
        warning('Adding the GeneratedBy fields failed...');
        fprintf('%s\n', ME.message);
    end

end


%% Add missing fields
function x = xASL_imp_FinishImport_AddMissingFields(x,dataPar)
    
    % Add fields that are in dataPar.x but missing in x
    if isfield(dataPar,'x')
        fieldsDataPar = fieldnames(dataPar.x);
        for iField = 1:numel(fieldsDataPar)
            if ~isfield(x,fieldsDataPar{iField,1}) && ~strcmp('dir',fieldsDataPar{iField,1})
                x.(fieldsDataPar{iField,1}) = dataPar.x.(fieldsDataPar{iField,1});
            end
        end
    end
    
end


