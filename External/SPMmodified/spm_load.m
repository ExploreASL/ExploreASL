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

% EXPLOREASL HACK: detect potentially erroneous empty cells
% Find all line ends
indexEol = find(S == eol);

if length(indexEol)<=1
	% In case the file does not have line ends (i.e. only a single line-end at the end of the file that is inserted automatically)
	% Then we consider that all data are provided on a single line. 
	% We then have to determine if the total number of cells can be divided by the number of columns
	
	% Read all cells
	d = textscan(S,'%s','Delimiter',delim);
	if rem(numel(d{1}),N) % can we divide cells/columns?
		error('Varying number of delimiters per line.');
	else
		d = reshape(d{1},N,[])'; % reshape
	end
else
	lineNumberEmpty = [];
	lineNumberIncorrectLength = [];
	
	% There are multiple line ends - we can parse per line
	% Create an empty cell array that will be filled line by line
	d = cell(length(indexEol), N);
	for iLine = 1:length(indexEol)
		% Read a single line from start to the next line end
		if iLine == 1
			lineD = textscan(S(1:indexEol(iLine)),'%s','Delimiter',delim);
		else
			lineD = textscan(S((indexEol(iLine-1)+1):indexEol(iLine)),'%s','Delimiter',delim);
		end
		
		% Count all cells and empty cells
		emptyCellCount = sum(cellfun(@(y) isempty(y), lineD{1}));
		cellCountPerLine = length(lineD{1});
		
		% Record empty lines
		if emptyCellCount > 0
			lineNumberEmpty = [lineNumberEmpty, iLine];
		end
		
		% Record lines of incorrect length
		if cellCountPerLine ~= N
			lineNumberIncorrectLength = [lineNumberIncorrectLength, iLine];
			
			% Read the cells from a single line to a table
			d(iLine, 1:min(cellCountPerLine,N)) = lineD{1}(1:min(cellCountPerLine,N));
			
			% Redefine type of empty cells from the initialized double to char
			for iFill = (min(cellCountPerLine,N)+1):N
				d{iLine, iFill} = '';
			end
		else
			d(iLine, :) = lineD{1};
		end
		
	end
end

%if ~isempty(lineNumberEmpty)
%    warning(['Found empty cells in ' f ' on lines ' num2str(lineNumberEmpty+1)]);
%end

if ~isempty(lineNumberIncorrectLength)
	warning(['Found lines with an incorrect length in ' f ' on lines ' num2str(lineNumberIncorrectLength)]);
	
	fprintf('%s\n', 'Repaired, but check carefully if this went OK');
	[fPath, fFile, fExt] = xASL_fileparts(f);
	pathSave = fullfile(fPath, [fFile '_repaired' fExt]);
	newCell = var;
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
