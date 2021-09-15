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



%% Compare field lists
function strError = compareFieldLists(jsonStructA,jsonStructB,fieldList)
    strError = '';
    
    % Threshold for the difference of numeric values
    threshNumeric = 1e-5;
    threshNumericArray = 1e-2;
    
    % Ignore version in Acknowledgements
    ignoreAcknowledgements = true;

    % Iterate over fields
    for iField=1:numel(fieldList)
        curFieldName = fieldList{iField};
        % Ignore version in Acknowledgements
		if ~ignoreAcknowledgements || ~strcmp(curFieldName,'Acknowledgements')
			% Otherwise continue comparing
			fieldContentA = jsonStructA.(fieldList{iField});
			fieldContentB = jsonStructB.(fieldList{iField});
			if isnumeric(fieldContentA) && isnumeric(fieldContentB)
				% Compare numbers
				if length(fieldContentA)==length(fieldContentB)
					if length(fieldContentA)==1
						% Compare numbers (check absolute difference)
						if abs(fieldContentA-fieldContentB)>threshNumeric
							strError = sprintf('%s           Different value: %s (%.6f vs %.6f)\n', strError,curFieldName,fieldContentA,fieldContentB);
						end
					else
						% Compare arrays (check sum of absolute differences)
						sumDiff = sum(abs(fieldContentA-fieldContentB));
						if sumDiff>threshNumericArray
							strError = sprintf('%s           Different value: %s (check arrays)\n', strError,curFieldName);
							% Set max number of elements to display
							maxNumElements = 10;
							if length(fieldContentA)<maxNumElements
								maxNumElements = length(fieldContentA);
							end
							% Initialize array view
							strError = sprintf('%s           [',strError);
							% Iterate over individual elements of array A
							for elField=1:length(fieldContentA)
								if elField<=maxNumElements
									if elField<maxNumElements
										strError = sprintf('%s%.4f, ',strError,fieldContentA(elField));
									elseif elField==length(fieldContentA)
										strError = sprintf('%s%.4f]\n',strError,fieldContentA(elField));
									elseif elField==maxNumElements
										strError = sprintf('%s%.4f ...]\n',strError,fieldContentA(elField));
									end
								end
							end
							% Initialize array view
							strError = sprintf('%s           [',strError);
							% Iterate over individual elements of array B
							for elField=1:length(fieldContentB)
								if elField<=maxNumElements
									if elField<maxNumElements
										strError = sprintf('%s%.4f, ',strError,fieldContentB(elField));
									elseif elField==length(fieldContentB)
										strError = sprintf('%s%.4f]\n',strError,fieldContentB(elField));
									elseif elField==maxNumElements
										strError = sprintf('%s%.4f ...]\n',strError,fieldContentB(elField));
									end
								end
							end
						end
					end
				else
					strError = sprintf('%s           Different dimension: %s\n', strError,curFieldName);
				end
			elseif ischar(fieldContentA) && ischar(fieldContentB)
				% Compare char arrays and strings
				if ~strcmp(fieldContentA,fieldContentB)
					strError = sprintf('%s           Different value: %s (%s vs %s)\n', strError,curFieldName,fieldContentA,fieldContentB);
				end
			elseif iscell(fieldContentA) && iscell(fieldContentB)
				% Compare cell arrays
				if ~(isempty(setdiff(fieldContentA,fieldContentB)) && isempty(setdiff(fieldContentB,fieldContentA)))
					strError = sprintf('%s           Different value: %s (array)\n', strError,curFieldName);
				end
			elseif isstruct(fieldContentA) && isstruct(fieldContentB)
				% Compare cell arrays
				if ~isequal(fieldContentA,fieldContentB)
					strError = sprintf('%s           Different value: %s (array)\n', strError,curFieldName);
				end
			elseif islogical(fieldContentA) && islogical(fieldContentB)
				% Compare numbers
				if ~(fieldContentA==fieldContentB)
					if fieldContentA
						fieldContentA = 'true';
					else
						fieldContentA = 'false';
					end
					
					if fieldContentB
						fieldContentB = 'true';
					else
						fieldContentB = 'false';
					end
					
					strError = sprintf('%s           Different value: %s (%s vs %s)\n', strError,curFieldName,fieldContentA,fieldContentB);
				end
			else
				% Neither number nor text
				if ~isequal(fieldContentA,fieldContentB)
					strError = sprintf('%s           Different value: %s (%s vs %s) - unknown or differing types\n', strError,curFieldName,fieldContentA,fieldContentB);
				end
			end
		end
    end
    
end





