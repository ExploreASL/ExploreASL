function xASL_adm_SaveJSON(data, jsonFileName)
% xASL_adm_SaveJSON Saves the values in the structure 'data' to a file in JSON format.
%
% Example:
%     data.name = 'chair';
%     data.color = 'pink';
%     data.metrics.height = 0.3;
%     data.metrics.width = 1.3;
%     xASL_adm_SaveJSON(data, 'out.json');
%
% Output 'out.json':
% {
% 	"name" : "chair",
% 	"color" : "pink",
% 	"metrics" : {
% 		"height" : 0.3,
% 		"width" : 1.3
% 		}
% 	}
%
%% WITH HACK HERE BY HM, for ExploreASL
% To allow cells, coming from BIDS JSON Decoding:
% QUICK & DIRTY FIX

    DataFields = fields(data);
    for iField=1:length(DataFields)
        CurrentField = data.(DataFields{iField});
        if iscell(CurrentField)
            data = rmfield(data,DataFields{iField});
            NewStruct = struct;
            NN=1;
            if isempty(CurrentField)
                NewStruct.(DataFields{iField}) = 'n/a';
            else
                for iX=1:size(CurrentField,1)
                    for iY=1:size(CurrentField,2)
                        for iZ=1:size(CurrentField,3)
                            NewStruct.([DataFields{iField} num2str(NN)]) = xASL_num2str(CurrentField{iX,iY,iZ});
                            NN=NN+1;
                        end
                    end
                end
            end
            data.(DataFields{iField}) = NewStruct;
        elseif isnumeric(CurrentField) && length(CurrentField)>1
            data = rmfield(data,DataFields{iField});
            NewStruct = struct;
            NN=1;
            for iX=1:size(CurrentField,1)
                for iY=1:size(CurrentField,2)
                    for iZ=1:size(CurrentField,3)
                        NewStruct.([DataFields{iField} num2str(NN)]) = CurrentField(iX,iY,iZ);
                        NN=NN+1;
                    end
                end
            end
            data.(DataFields{iField}) = NewStruct;
        end
    end

    fclose all;
    xASL_adm_CreateDir(fileparts(jsonFileName));
    fid = fopen(jsonFileName,'w');
    %fid=1;
    for iD=1:length(data)
        writeElement(fid, data(iD),'');
    end
    fprintf(fid,'\n');
    fclose(fid);
end

function writeElement(fid, data,tabs)
    namesOfFields = fieldnames(data);
    numFields = length(namesOfFields);
    tabs = sprintf('%s\t',tabs);
    fprintf(fid,'{\n%s',tabs);
   

    for i = 1:numFields - 1
        currentField = namesOfFields{i};
        currentElementValue = data.(currentField);
        
        if isnumeric(currentElementValue) && length(currentElementValue)>1 %% ExploreASL HACK
            currentElementValue = xASL_num2str(currentElementValue);
        end
        
        writeSingleElement(fid, currentField,currentElementValue,tabs);
        fprintf(fid,',\n%s',tabs);
    end
    if isempty(i)
        i=1;
    else
      i=i+1;
    end
      
    
    currentField = namesOfFields{i};
    currentElementValue = data.(currentField);
    writeSingleElement(fid, currentField,currentElementValue,tabs);
    fprintf(fid,'\n%s}',tabs);
end

function writeSingleElement(fid, currentField,currentElementValue,tabs)
    
        % if this is an array and not a string then iterate on every
        % element, if this is a single element write it
        if length(currentElementValue) > 1 && ~ischar(currentElementValue)
            fprintf(fid,' "%s" : [\n%s',currentField,tabs);
            for m = 1:length(currentElementValue)-1
                writeElement(fid, currentElementValue(m),tabs);
                fprintf(fid,',\n%s',tabs);
            end
            if isempty(m)
                m=1;
            else
              m=m+1;
            end
            
            writeElement(fid, currentElementValue(m),tabs);
          
            fprintf(fid,'\n%s]\n%s',tabs,tabs);
        elseif isstruct(currentElementValue)
            fprintf(fid,'"%s" : ',currentField);
            writeElement(fid, currentElementValue,tabs);
        elseif isnumeric(currentElementValue)
            fprintf(fid,'"%s" : %g' , currentField,currentElementValue);
        elseif isempty(currentElementValue)
            %fprintf(fid,'"%s" : null' , currentField,currentElementValue);
			fprintf(fid,'"%s" : ""', currentField);
        else %ischar or something else ...
            fprintf(fid,'"%s" : "%s"' , currentField,currentElementValue);
        end
end