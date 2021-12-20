function strError = xASL_bids_CompareFieldLists(jsonStructA, jsonStructB, fieldList, ignoreFields)
%xASL_bids_CompareFieldLists Compare JSON file field lists
%
% FORMAT: strError = xASL_bids_CompareFieldLists(jsonStructA, jsonStructB, fieldList, ignoreFields)
%
% INPUT:
%         jsonStructA - struct of JSON A (STRUCT, REQUIRED)
%         jsonStructB - struct of JSON A (STRUCT, REQUIRED)
%         fieldList   - List of shared fields (LIST, REQUIRED)
%         ignoreFields - List of fields which should be ignored (OPTIONAL, DEFAULT = none)
%
% OUTPUT:
%         strError    - Error string
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      This script compares the content of two JSON files for
%                   the BIDS flavor testing.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          ...
%
% REFERENCES:       ...
% __________________________________
% Copyright (c) 2015-2021 ExploreASL

    if nargin < 4
        ignoreFields = {};
    end

    strError = '';
    
    % Threshold for the difference of numeric values
    threshNumeric = 1e-5;
    threshNumericArray = 1e-2;

    % Iterate over fields
    for iField=1:numel(fieldList)
        doCompare = true;
        curFieldName = fieldList{iField};
        fieldContentA = jsonStructA.(fieldList{iField});
        fieldContentB = jsonStructB.(fieldList{iField});
        % Check if field is supposed to be ignored
        if ~isempty(find(ismember(ignoreFields, curFieldName), 1))
            doCompare = false;
        end
        if doCompare
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


