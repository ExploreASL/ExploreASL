function [IsValid] = xASL_qc_CheckValidityJSON( PathJSON )
%xASL_qc_CheckValidityJSON This function loads a QC JSON (simply JSON file, won't
%take any exotic files) and simply check whether there is any empty value
%after a key. If this is the case, it will throw a warning, which will skip
%reading this JSON by the compiled spm_jsonread, avoiding the crash that
%this may result in.

IsValid = true;

CellData = importdata(PathJSON);

% we loop across all cells, trying to recreate the JSON in a struct

% First we search for the first layer
% IndicesFirst = cellfun(@(y) strcmp(y,'{'), CellData));
% IndicesFirst2 = cellfun(@(y) strcmp(y,'}'), CellData);

for iCell=1:size(CellData,1)
    
    CurrentCell = CellData{iCell,1};
        
    IndicesQuotes = strfind(CurrentCell,'"');
    IndexColon = strfind(CurrentCell,':');
%     
%     Field1 = []; % we start with no first field
%     Field2 = []; % we start with no second field
    
    if ~isempty(IndicesQuotes) && ~isempty(IndexColon)
        % only now continue, we need a key & value
        % if IndicesQuotes(1)==2 % first field
        % if IndicesQuotes(1)==3 % second field -> etc, we could use this
        % to "repair" the JSON. But for now we use this script only to
        % quickly check the JSON validity
        
        Key = CurrentCell(IndicesQuotes(1)+1:IndicesQuotes(2)-1);
        Value = CurrentCell(IndexColon(1)+1:end);
        if strcmp(Value(end),',')
            Value = Value(1:end-1);
        end
        
        if isempty(Value) || strcmp(Value,' ') % find empty values
            if IsValid % if we didnt throw the warning already
                warning(['Invalid JSON: ' PathJSON]);
            end
            IsValid = false;
%             Value = NaN;
            

%         else
%             Value = Value(1:end-1);
        end
    end

end

