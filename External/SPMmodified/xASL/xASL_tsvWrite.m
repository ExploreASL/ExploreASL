function xASL_tsvWrite(InputCell, PathTSV, bOverwrite, bCSV)
%xASL_tsvWrite Write cell array to TSV
%
% FORMAT: xASL_tsvWrite(InputCell, PathTSV[, bOverwrite, bCSV])
%
% INPUT:
%   InputCell           - cell array, containing a combination of numerical
%                         values/strings (REQUIRED)
%   PathTSV             - output file, with TSV extension (REQUIRED)
%   bOverwrite          - boolean to specify if an existing file should be
%                         overwritten (OPTIONAL, DEFAULT=false)
%   bCSV                - boolean stating if we create a CSV (1) instead of
%                         TSV file (0) (OPTIONAL, DEFAULT=false)
% OUTPUT: n/a 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function loads a cell array and prints it to a
% tab-separated value (TSV) file, which is the format that BIDS prefers.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_tsvWrite(ParticipantsMetadata, '/MyStudy/participants.tsv');
% __________________________________
% Copyright 2015-2020 ExploreASL

%% -------------------------------------------------------
%% Admin
if nargin<4 || isempty(bCSV) || bCSV==0
	bCSV = 0;
    DelimiterIs = '%s\t'; % tab-separated values (TSV)
elseif bCSV==1
    DelimiterIs = '%s,'; % comma-separated values (CSV)
else
    warning('Wrong choice of delimiter, bCSV should be CSV (1) or TSV (0)');
    return;
end
if nargin<3 || isempty(bOverwrite)
    bOverwrite = false;
end
if nargin<2 || isempty(InputCell) || isempty(PathTSV)
    error('Specify both InputCell and PathTSV');
end

[Fpath, Ffile, Fext] = xASL_fileparts(PathTSV);
if ~strcmp(Fext,'.tsv') && ~strcmp(Fext,'.csv')
    warning('Wrong extension, resetting to default (.tsv).');
    PathTSV = fullfile(Fpath, [Ffile '.tsv']);
end

% Extension was set to .tsv but bCSV is true or the other way round
if (bCSV==1 && strcmp(Fext,'.tsv')) || (bCSV==0 && strcmp(Fext,'.csv'))
    error('File extension and selection of tsv/csv do not match...');
end

if xASL_exist(PathTSV, 'file') && ~bOverwrite
    warning([PathTSV ' already existed, skipping']);
    return;
elseif bOverwrite
    xASL_delete(PathTSV);
end


%% -------------------------------------------------------
%% Create the file
xASL_adm_CreateDir(fileparts(PathTSV)); % avoid errors with non-existing folder
FileID = fopen(PathTSV ,'w+'); % create a new file for writing
if FileID<0 % give warning when cannot create the file
    fprintf(['Filename is ' PathTSV '\n']);
    warning('Is something wrong with the path we try to save the file, perhaps it is too long?');
end


%% Iterate through cells with iX and iY coordinates & print the cell array contents
for iX=1:size(InputCell,1)
    for iY=1:size(InputCell,2)
        if iY==size(InputCell,2)
            fprintf(FileID,'%s', xASL_num2str(InputCell{iX,iY}));
        else
            fprintf(FileID,DelimiterIs, xASL_num2str(InputCell{iX,iY}));
        end
    end
    if iX<size(InputCell,1)
        fprintf(FileID,'\n');
    end
end


%% Close the file & unlock it
fclose(FileID); 


end
