function [CellContents] = xASL_tsvRead(PathTSV)
%xASL_tsvRead Read a TSV file and output contents to cell array
%
% FORMAT: [CellContents] = xASL_tsvRead(PathTSV)
%
% INPUT:
%   PathTSV             - file, with TSV extension (REQUIRED)
%
% OUTPUT:
%   CellContents        - cell array, containing a combination of numerical
%                         values/strings
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function loads a tab-separated value (TSV) file - which
% is the format that BIDS prefers - and outputs it to a cell array.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: ParticipantsMetadata = xASL_tsvRead('/MyStudy/participants.tsv');
% __________________________________
% Copyright 2015-2020 ExploreASL


%% -------------------------------------------------------
%% Admin
if nargin<1 || isempty(PathTSV)
    warning('Invalid PathTSV input argument, skipping');
    return;
end    

[Fpath, Ffile, Fext] = xASL_fileparts(PathTSV);
if ~strcmp(Fext,'.tsv')
    warning('Wrong extension, resetting to .tsv');
    PathTSV = fullfile(Fpath, [Ffile '.tsv']);
end

%% -------------------------------------------------------
%% Read the file
ReadCell = textread(PathTSV, '%s', 'delimiter','\t\n', 'bufsize', 10000000);
% now we have a 1D cell structure, where the \n became empty cells (other stuff is skipped)
iRow = 1;
iColumn = 1;
for iCell=1:length(ReadCell)
    if isempty(ReadCell{iCell})
        iRow = iRow+1;
        iColumn = 1;
    else
        CellContents{iRow,iColumn} = ReadCell{iCell};
        iColumn = iColumn+1;
    end
end

    
    

end

