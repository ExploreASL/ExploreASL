function x = spm_load(f,v)
% Load text and numeric data from file
% FORMAT x = spm_load(f,v)
% f  - filename (can be gzipped) {txt,mat,csv,tsv,json}
% v  - name of field to return if data stored in a structure [default: '']
%      or index of column if data stored as an array
%
% x  - corresponding data array or structure
%__________________________________________________________________________
% Copyright (C) 1995-2017 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_load.m 7097 2017-06-07 13:53:55Z guillaume $


%-Get a filename if none was passed
%--------------------------------------------------------------------------
if ~nargin
    [f,sts] = spm_select(1,{...
        'mat',...                        % *.txt, *.mat
        '^.*\.csv$','^.*\.csv.gz$',...   % *.csv, *.csv.gz
        '^.*\.tsv$','^.*\.tsv.gz$',...   % *.tsv, *.tsv.gz
        '^.*\.json$','^.*\.json.gz$',... % *.json, *.json.gz
        });
    if ~sts, x = []; return; end
end

if ~exist(f,'file')
    error('Unable to read file ''%s''',f);
end

if nargin < 2, v = ''; end

%-Load the data file
%--------------------------------------------------------------------------
switch spm_file(f,'ext')
    case 'txt'
        x = load(f,'-ascii');
    case 'mat'
        x  = load(f,'-mat');
    case 'csv'
        % x = csvread(f); % numeric data only
        x = dsvread(f,',');
    case 'tsv'
        % x = dlmread(f,'\t'); % numeric data only
        x = dsvread(f,'\t');
    case 'json'
        x = spm_jsonread(f);
    case 'gz'
        fz  = gunzip(f,tempname);
        sts = true;
        try
            x   = spm_load(fz{1});
        catch
            sts = false;
        end
        delete(fz{1});
        rmdir(spm_file(fz{1},'path'));
        if ~sts, error('Cannot load ''%s''.',f); end
    otherwise
        try
            x = load(f);
        catch
            error('Unknown file format.');
        end
end

%-Return relevant subset of the data if required
%--------------------------------------------------------------------------
if isstruct(x)
    if isempty(v)
        fn = fieldnames(x);
        if numel(fn) == 1 && isnumeric(x.(fn{1}))
            x = x.(fn{1});
        end
    else
        if ischar(v)
            try
                x = x.(v);
            catch
                error('Data do not contain array ''%s''.',v);
            end
        else
            fn = fieldnames(x);
            try
                x = x.(fn{v});
            catch
                error('Invalid data index.');
            end
        end
    end
elseif isnumeric(x)
    if isnumeric(v)
        try
            x = x(:,v);
        catch
            error('Invalid data index.');
        end
    elseif ~isempty(v)
        error('Invalid data index.');
    end
end


end % EXPLORE HACK

%==========================================================================
% function x = dsvread(f,delim)
%==========================================================================
function x = dsvread(f,delim)
% Read delimiter-separated values file into a structure array
%  * header line of column names will be used if detected
%  * 'n/a' fields are replaced with NaN

%-Input arguments
%--------------------------------------------------------------------------
if nargin < 2, delim = '\t'; end
delim = sprintf(delim);
eol   = sprintf('\n');

%-Read file
%--------------------------------------------------------------------------
S   = fileread(f);
if isempty(S), x = []; return; end
if S(end) ~= eol, S = [S eol]; end
S   = regexprep(S,{'\r\n','\r','(\n)\1+'},{'\n','\n','$1'});

%-Get column names from header line (non-numeric first line)
%--------------------------------------------------------------------------
h   = find(S == eol,1);
hdr = S(1:h-1);
var = regexp(hdr,delim,'split');

% EXPLOREASL HACK: manage trailing \t on header only
if isempty(var{end})
    var = var(1:end-1);
end

N   = numel(var);
n1  = isnan(cellfun(@str2double,var));
n2  = cellfun(@(x) strcmpi(x,'NaN'),var);
if any(n1 & ~n2)
    hdr     = true;
    try
        var = genvarname(var);
    catch
        var = matlab.lang.makeValidName(var,'ReplacementStyle','hex');
        var = matlab.lang.makeUniqueStrings(var);
    end
    S       = S(h+1:end);
else
    hdr     = false;
    var     = cell(N,1);
    for i=1:N
        var{i} = sprintf(['var%0' num2str(floor(log10(N))+1) 'd'],i);
    end
end

%-Parse file
%--------------------------------------------------------------------------
if strcmpi(spm_check_version,'octave') % bug #51093
    S = strrep(S,delim,'#');
    delim = '#';
end
d = textscan(S,'%s','Delimiter',delim);

% EXPLOREASL HACK: detect potentially erroneous empty cells
illegalEmptyCells = DetectIllegalEmptyCells(d, N);

if ~isempty(illegalEmptyCells)
    warning(['Found illegal empty cells in ' f ', trying to repair']);
    d = FixIllegalEmptyCells(d, N);

    % Check if it was repaired
    if rem(numel(d{1}),N) % can we divide cells/columns?
        fprintf('%s\n', ['Could not repair ' f]);
    else 
        % try to save it repaired file
        fprintf('%s\n', ['Could repair ' f ', but check carefully if this went OK']);
        [fPath, fFile, fExt] = xASL_fileparts(f);
        pathSave = fullfile(fPath, [fFile '_repaired' fExt]);
        newCell = var;
        d = reshape(d{1},N,[])'; % reshape
        newCell(2:size(d,1)+1,:) = d;

        if strcmp(fExt, '.tsv')
            fprintf('%s\n', ['Saved repaired file: ' pathSave]);
            xASL_tsvWrite(newCell, pathSave);
        elseif strcmp(fExt, '.csv')
            fprintf('%s\n', ['Saving repaired file: ' pathSave]);
            xASL_csvWrite(newCell, pathSave);
        else
            fprintf('%s\n', ['Unknown delimiter (extension: ' fExt ', could not save repaired file']);
        end
    end
elseif rem(numel(d{1}),N) % can we divide cells/columns?
    error('Varying number of delimiters per line.');
else
    d = reshape(d{1},N,[])'; % reshape
end



allnum = true;
for i=1:numel(var)
    sts = true;
    dd = zeros(size(d,1),1);
    for j=1:size(d,1)
        if strcmp(d{j,i},'n/a')
            dd(j) = NaN;
        else
            dd(j) = str2double(d{j,i}); % i,j considered as complex
            if isnan(dd(j)), sts = false; break; end
        end
    end
    if sts
        x.(var{i}) = dd;
    else
        x.(var{i}) = d(:,i);
        allnum     = false;
    end
end

if ~hdr && allnum
    x = struct2cell(x);
    x = [x{:}];
end


end

%% ========================================================================

function [d] = FixIllegalEmptyCells(d, N)
%FixIllegalEmptyCells % EXPLOREASL HACK: Fix potentially erroneous empty cells
%   d = all cells from input delimiter-separate value file (without header)
%   N = number of columns
    
    illegalEmptyCells = DetectIllegalEmptyCells(d, N);

    while ~isempty(illegalEmptyCells)
        illegalRow = illegalEmptyCells(1);
        fprintf('%s\n', ['Found empty cells at row ' xASL_num2str(illegalRow)]);

        % Remove all adjacent empty cells
        while isempty(d{1}{illegalRow})
            d{1}(illegalRow:end-1) = d{1}(illegalRow+1:end);
            d{1} = d{1}(1:end-1);
        end

        illegalEmptyCells = DetectIllegalEmptyCells(d, N);
    end

end

%% ========================================================================


function [illegalEmptyCells] = DetectIllegalEmptyCells(d, N)
%DetectIllegalEmptyCells % EXPLOREASL HACK: detect & fix potentially erroneous empty cells
%   d = all cells from input delimiter-separate value file (without header)
%   N = number of columns
    
    emptyCells = find(cellfun(@(y) isempty(y), d{1}));
    % N = number of columns
    % so if any index-1 of the empty cells equals to the number of columns,
    % this means that the empty cell is at the end of a line and potentially
    % crashing this function trying to open a delimiter-separated value file
    % (dsv)
    
    illegalEmptyCells = (emptyCells-1)/N == round((emptyCells-1)/N);
    illegalEmptyCells = emptyCells(illegalEmptyCells);

end