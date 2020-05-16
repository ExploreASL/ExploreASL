function xASL_tsvWrite(InputCell, PathTSV, bOverwrite)
%xASL_tsvWrite Write cell array to TSV
%
% FORMAT: xASL_adm_tsvWrite(InputCell, PathTSV[, bOverwrite])
%
% INPUT:
%   InputCell           - cell array, containing a combination of numerical
%                         values/strings (REQUIRED)
%   PathTSV             - output file, with TSV extension (REQUIRED)
%   bOverwrite          - boolean to specify if an existing file should be
%                         overwritten (OPTIONAL, DEFAULT=false)
% OUTPUT: n/a 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function loads a cell array and prints it to a
% tab-separated value (TSV) file, which is the format that BIDS prefers.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_adm_tsvWrite(ParticipantsMetadata, '/MyStudy/participants.tsv');
% __________________________________
% Copyright 2015-2020 ExploreASL

%% -------------------------------------------------------
%% Admin

if nargin<3 || isempty(bOverwrite)
    bOverwrite = false;
end
if nargin<2 || isempty(InputCell) || isempty(PathTSV)
    error('Specify both InputCell and PathTSV');
end

[Fpath, Ffile, Fext] = xASL_fileparts(PathTSV);
if ~strcmp(Fext,'.tsv')
    warning('Wrong extension, resetting to .tsv');
    PathTSV = fullfile(Fpath, [Ffile '.tsv']);
end

if xASL_exist(PathTSV, 'file') && ~bOverwrite
    warning([PathTSV ' already existed, skipping']);
    return;
elseif bOverwrite
    xASL_delete(PathTSV);
end

%% -------------------------------------------------------
%% Create the file

FileID = fopen(PathTSV ,'w+'); % create a new file for writing

% iterate through cells with iX and iY coordinates
% & print the cell array contents
for iX=1:size(InputCell,1)
    for iY=1:size(InputCell,2)
        fprintf(FileID,'%s\t', xASL_num2str(InputCell{iX,iY}));
    end
    fprintf(FileID,'\n');
end

fclose(FileID); % stop writing to the file, close & unlock it 


end


