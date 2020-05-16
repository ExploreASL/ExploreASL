function [CellContents] = xASL_csvRead(PathCSV)
%xASL_csvRead Read a csv file and output contents to cell array
%
% FORMAT: [CellContents] = xASL_csvRead(PathCSV)
%
% INPUT:
%   PathCSV             - file, with csv extension (REQUIRED)
%
% OUTPUT:
%   CellContents        - cell array, containing a combination of numerical
%                         values/strings
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function loads a comma-separated value (csv) file - which
% is the format that BIDS prefers - and outputs it to a cell array.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: ParticipantsMetadata = xASL_csvRead('/MyStudy/participants.csv');
% __________________________________
% Copyright 2015-2020 ExploreASL


%% -------------------------------------------------------
%% Admin
if nargin<1 || isempty(PathCSV)
    warning('Invalid Pathcsv input argument, skipping');
    return;
end    

[Fpath, Ffile, Fext] = xASL_fileparts(PathCSV);
if ~strcmp(Fext,'.csv')
    warning('Wrong extension, resetting to .csv');
    PathCSV = fullfile(Fpath, [Ffile '.csv']);
end

%% -------------------------------------------------------
%% Read the file

TempText = textread(PathCSV, '%s', 'delimiter','\n', 'bufsize', 10000000);

for iT=1:size(TempText,1)
    CommaIndices = find(TempText{iT,1}==',');
    CommaIndices = [0 CommaIndices length(TempText{iT,1})+1];
    for iC=1:length(CommaIndices)-1
        CI1 = CommaIndices(iC)+1;
        CI2 = CommaIndices(iC+1)-1;
        CellContents{iT,iC} = TempText{iT,1}(CI1:CI2);

        %% Replaces spaces by underscores
        iSpace              = find(CellContents{iT,iC}==' ');
        while ~isempty(iSpace)
            if iSpace(1)==1
                   CellContents{iT,iC} = CellContents{iT,iC}(2:end);
            elseif iSpace(1)==length(CellContents{iT,iC})
                   CellContents{iT,iC} = CellContents{iT,iC}(1:end-1);
            else
                CellContents{iT,iC} = [CellContents{iT,iC}(1:iSpace(1)-1) '_' CellContents{iT,iC}(iSpace(1)+1:end)];
            end
            iSpace = find(CellContents{iT,iC}==' ');
        end

    end
end
    
    

end

