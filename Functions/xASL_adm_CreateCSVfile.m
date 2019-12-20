function xASL_adm_CreateCSVfile(CSVfilename,CSVdata)
%xASL_adm_CreateCSVfile Creates a CSV file that can be opened with excel from
%your data
% CSVfilename=FilePathName
% CSVdata = your data

    fileID      = fopen(CSVfilename,'w');
    for iY=1:size(CSVdata,1)
        for iX=1:size(CSVdata,2)
            fprintf(fileID,CSVdata{iY,iX});
            if  iX<size(CSVdata,2); fprintf(fileID,','); end
        end
        if  iY<size(CSVdata,1); fprintf(fileID,'\n'); end
    end
    fclose(fileID);

end



        
