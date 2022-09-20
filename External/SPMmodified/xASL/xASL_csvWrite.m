function xASL_csvWrite(InputCell, PathCSV, bOverwrite)
%xASL_csvWrite Write cell array to CSV
%
% FORMAT: xASL_csvWrite(InputCell, PathCSV, bOverwrite)
%
% INPUT:    InputCell   - cell array
%           PathCSV     - path to file
%           bOverwrite  - overwrite file
%
% OUTPUT:   n/a
%
%
% DESCRIPTION: Rudimentary function, please use xASL_tsvWrite instead.
% For usage, type help xASL_tsvWrite.
% This function will still work though.
%
% EXAMPLE: n/a
%
% __________________________________
% Copyright 2015-2020 ExploreASL


%% -------------------------------------------------------
%% Admin
if nargin<3 || isempty(bOverwrite)
    bOverwrite = false;
end
if nargin<2 || isempty(InputCell) || isempty(PathCSV)
    error('Specify both InputCell and PathTSV');
end

warning('You are creating a CSV file, please consider using the TSV-file format, which is default in ExploreASL (per BIDS)');
bCSV = 1;
xASL_tsvWrite(InputCell, PathCSV, bOverwrite, bCSV);


end



        
