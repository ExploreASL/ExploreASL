function [CellContents] = xASL_tsvRead(PathTSV, bStruct)
%xASL_tsvRead Read a TSV file and output contents to cell array
%
% FORMAT: [CellContents] = xASL_tsvRead(PathTSV)
%
% INPUT:
%   PathTSV             - file, with TSV extension (REQUIRED)
%   bStruct             - parse to a struct with separate field per tsv
%                         column (OPTIONAL, DEFAULT = false)
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
if nargin<2 || isempty(bStruct)
    bStruct = false;
end

[Fpath, Ffile, Fext] = xASL_fileparts(PathTSV);
if ~strcmp(Fext,'.tsv')
    warning('Wrong extension, resetting to .tsv');
    PathTSV = fullfile(Fpath, [Ffile '.tsv']);
end

%% -------------------------------------------------------
%% Read the file into multiple cell structures
ReadCell = spm_load(PathTSV);

if bStruct
    % Keep output as structs
    CellContents = ReadCell;
else
    % Parse struct and convert to cell array
    FieldsAre = fields(ReadCell);
    for iField=1:length(FieldsAre)
        % Header
        CellContents{1, iField} = FieldsAre{iField};

        % Content
        CurrentContents = ReadCell.(FieldsAre{iField});

        if ~iscell(CurrentContents)
            for iRow=1:length(CurrentContents)
                NewContents{iRow, 1} = CurrentContents(iRow, 1);
            end
        else
            NewContents = CurrentContents;
        end

        CellContents(2:length(NewContents)+1, iField) = NewContents;
    end
end
    
    
end