function xASL_csvWrite(InputCell, PathCSV, bOverwrite)
%xASL_csvWrite Write cell array to CSV
%
% FORMAT: xASL_csvWrite(InputCell, PathCSV[, bOverwrite])
%
% INPUT:    InputCell   - cell array (REQUIRED)
%           PathCSV     - path to file (REQUIRED)
%           bOverwrite  - overwrite file (OPTIONAL, DEFAULT = false)
%
% OUTPUT:   n/a
%
%
% DESCRIPTION: Legacy function, please use xASL_tsvWrite instead, as this function might be discontinued in the future.
% For usage, type help xASL_tsvWrite.
%
% EXAMPLE: n/a
% __________________________________
% Copyright 2015-2022 ExploreASL

%% -------------------------------------------------------
%% Admin
if nargin<3 || isempty(bOverwrite)
    bOverwrite = false;
end
if nargin<2 || isempty(InputCell) || isempty(PathCSV)
    error('Specify both InputCell and PathCSV');
end

warning('You are creating a CSV file, please consider using the TSV-file format, which is default in ExploreASL (per BIDS)');
bCSV = 1;
xASL_tsvWrite(InputCell, PathCSV, bOverwrite, bCSV);

end
