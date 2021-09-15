function [differences,identical,dn] = xASL_bids_CompareStructuresJSON(differences,identical,bPrintReport,allFiles,iFile,dn,currentFileA,currentFileB)
%xASL_bids_CompareStructuresJSON Compare JSON files
%
% FORMAT: [differences,identical,dn] = xASL_bids_CompareStructuresJSON(differences,identical,bPrintReport,allFiles,iFile,dn,currentFileA,currentFileB)
%
% INPUT:
%         ...
%
% OUTPUT:
%         ...
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      ...
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          ...
%
% REFERENCES:       ...
% __________________________________
% Copyright @ 2015-2021 ExploreASL


    % Ignore Acknowledgements (setting to ignore version differences (ExploreASL v1.x.x) there)
    ignoreAcknowledgements = true;

    jsonErrorReport='';

    % Compare JSON files on field basis
    if (exist(currentFileA,'file') && exist(currentFileB,'file')) % xASL_exist somehow didn't work here (again)
        % Import JSON files
        jsonA = spm_jsonread(char(currentFileA));
        jsonB = spm_jsonread(char(currentFileB));
        
        % Make all paths to paths of THIS os
        jsonA = fixPathFields(jsonA);
        jsonB = fixPathFields(jsonB);

        % Get JSON field names
        fieldNamesA = fieldnames(jsonA);
        fieldNamesB = fieldnames(jsonB);

        % Check which fields are shared and which different
        sharedFieldsAB = intersect(fieldNamesB,fieldNamesA);

        % Fields that are in B, but missing in A
        missingFields = setdiff(fieldNamesB,fieldNamesA);
        % Print out missing fields
        for iField=1:numel(missingFields)
            if ignoreAcknowledgements && strcmp(missingFields{iField},'Acknowledgements')
                continue
            end
            jsonErrorReport = sprintf('%s           Missing field: %s\n',jsonErrorReport,missingFields{iField});
        end

        extraFields = setdiff(fieldNamesA,fieldNamesB);
        % Print out missing fields
        for iField=1:numel(extraFields)
            if ignoreAcknowledgements && strcmp(extraFields{iField},'Acknowledgements')
                continue
            end
            jsonErrorReport = sprintf('%s           Extra field: %s\n',jsonErrorReport,extraFields{iField});
        end

        % Now we can compare these fields like in the part above
        jsonErrorReport = [jsonErrorReport, compareFieldLists(jsonA,jsonB,sharedFieldsAB)];

        if ~isempty(jsonErrorReport)
            if bPrintReport
                fprintf('File:      %s\n',allFiles{iFile});
                fprintf('%s',jsonErrorReport);
            end

            % Save difference
            differences{dn,1} = ['Different file content: ', allFiles{iFile}, ' '];
            dn = dn+1;
        end

    end

end



%% Fix path fields
function jsonStruct = fixPathFields(jsonStruct)

    fieldNames = fieldnames(jsonStruct);
    numOfField = size(fieldNames,1);
    
    % Check number of fields
    if numOfField > 0
        for iField = 1:numOfField
            currentValue = jsonStruct.(fieldNames{iField,1});
            if ischar(currentValue)
                % Convert all paths to unix paths
                jsonStruct.(fieldNames{iField,1}) = strrep(currentValue,'\','/');
            end
        end
    end

end




