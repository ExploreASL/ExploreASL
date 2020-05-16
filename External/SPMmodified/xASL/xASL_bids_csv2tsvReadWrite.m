function [PathTSV, CellContents] = xASL_bids_csv2tsvReadWrite(PathIn, bDeleteCSV, bWriteTSV)
%xASL_bids_csv2tsvReadWrite Multi-purpose read & convert csv/tsv
%
% FORMAT: [PathTSV, CellContents] = xASL_bids_csv2tsvReadWrite(PathIn[, bDeleteCSV, bWriteTSV])
%
% INPUT:
%   PathIn          - path to CSV or TSV input file (REQUIRED)
%   bDeleteCSV      - boolean specifying if an original CSV should be deleted (OPTIONAL, DEFAULT=true)
%   bOverwrite      - boolean specifying if the CSV contents should be written to a new TSV file (OPTIONAL, DEFAULT=true)
%
% OUTPUT:
%   PathTSV         - path to TSV file
%   CellContents    - cell array, containing a combination of numerical values/strings
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function PathIn and loads it, also trying CSV or TSV
% extensions if these exist. It outputs the contents to a cell array. If a
% CSV file exists but not a TSV file, it converts and replaces the CSV to
% TSV file, per BIDS. This function has the following parts:
%
% 1) Read the CSV or TSV file
% 2) Write the TSV file (if requested)
% 3) Delete the CSV file (if requested)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_adm_tsvWrite(ParticipantsMetadata, '/MyStudy/participants.tsv');
% __________________________________
% Copyright 2015-2020 ExploreASL



%% -------------------------------------------------------
%% Admin
if nargin<3 || isempty(bWiteTSV)
    bWriteTSV = true;
end

if nargin<2 || isempty(bDeleteCSV)
    bDeleteCSV = true;
end

if nargin<1 || isempty(PathIn)
    warning('Invalid input path argument, skipping');
    return;
end

[Fpath, Ffile, Fext] = xASL_fileparts(PathIn);

if isempty(regexp(Fext, '^\.(tsv|csv)$'))
    warning('Wrong extension, resetting to .csv');
    PathIn = fullfile(Fpath, [Ffile '.csv']);
end

PathTSV = fullfile(Fpath,[Ffile '.tsv']);

if strcmp(Fext, '.csv') && ~exist(PathIn, 'file')
    % if we want to read a CSV file, but it doesnt exist, try reading a TSV file instead
    Fext = '.tsv';
    if exist(PathTSV, 'file')
        PathIn = PathTSV;
    else
        warning([PathIn ' didnt exist, skipping']);
        return;
    end
end


%% -------------------------------------------------------
%% 1) Read the CSV or TSV file
if strcmp(Fext,'.tsv')
    CellContents = xASL_tsvRead(PathIn);
elseif strcmp(Fext,'.csv')
    CellContents = xASL_csvRead(PathIn);
end


%% -------------------------------------------------------
%% 2) Write the TSV file (if requested)
if ~exist(PathTSV, 'file') && bWriteTSV
    xASL_tsvWrite(CellContents, PathTSV, 1); % overwrite
end


%% -------------------------------------------------------
%% 3) Delete the CSV file (if requested)
if bDeleteCSV && exist(PathTSV, 'file') && strcmp(PathIn(end-2:end), 'csv')
    xASL_delete(PathIn);
end
   

end