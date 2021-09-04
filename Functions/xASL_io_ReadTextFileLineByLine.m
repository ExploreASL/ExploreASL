function textArray = xASL_io_ReadTextFileLineByLine(pathTextFile)
%xASL_io_ReadTextFileLineByLine Read a text file line by line
%
% FORMAT: textArray = xASL_io_ReadTextFileLineByLine(pathTextFile)
% 
% INPUT:
%   pathTextFile   - Path to file  (CHAR ARRAY, REQUIRED)
%
% OUTPUT:
%   textArray      - Cell array containing file text
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Read a text file line by line.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        textArray = xASL_io_ReadTextFileLineByLine(pathTextFile);
% __________________________________
% Copyright 2015-2021 ExploreASL

    % Open text file
    fid = fopen(pathTextFile);
    textLine = fgetl(fid);
    iLine = 1;
    
    % Read text file line by line
    while ischar(textLine)
        textArray{iLine,1} = textLine;
        textLine = fgetl(fid);
        iLine = iLine+1;
    end
    
    % Close text file
    fclose(fid);


end


