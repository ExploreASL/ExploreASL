function xASL_io_MergeTextFiles(pathA,pathB,pathOut,headerText)
%xASL_io_MergeTextFiles Merge text files A and B and write the output to the pathOut file
%
% FORMAT: xASL_io_MergeTextFiles(pathA,pathB,pathOut[,headerText])
% 
% INPUT:
%   pathA        - Path to file A (CHAR ARRAY, REQUIRED)
%   pathB        - Path to file B (CHAR ARRAY, REQUIRED)
%   pathOut      - Path to output file (CHAR ARRAY, REQUIRED)
%   headerText   - Header text (CHAR ARRAY, OPTIONAL)
%
% OUTPUT:
%   n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Merge text files A and B and write the output to the pathOut file.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        xASL_io_MergeTextFiles(pathA,pathB,pathOut,headerText);
% __________________________________
% Copyright 2015-2021 ExploreASL

    % Input check
    if nargin < 4 || isempty(headerText)
        headerText = '';
    end

    % Read text files
    textA = xASL_io_ReadTextFileLineByLine(pathA);
    textB = xASL_io_ReadTextFileLineByLine(pathB);
    
    % Write text to file
    fid = fopen(pathOut,'wt');
    if ~isempty(headerText)
        fprintf(fid, '%s\n', headerText);
    end
    for iLine = 1:size(textA,1)
        fprintf(fid, '%s\n', textA{iLine,1});
    end
    for iLine = 1:size(textB,1)
        fprintf(fid, '%s\n', textB{iLine,1});
    end
    fclose(fid);


end



