function [tree, optionalTokens] = xASL_adm_FindByRegExp(root, dirSpecs, varargin)
% FindByPath find files or folders recursively by specifying search patterns
%
% FORMAT: 
%     xasl_adm_FindByRegExp(root, dirSpecs[, varargin])
% INPUT:
%     root      - root directory to start search
%     dirSpecs  - string, cellstr or structure containing regular expressions for recursive search
%                 Examples:
%                 - simple string:
%                     '\d{3}|^PCASL.*|PCASL.*\.PAR$'           (i.e., use '|' as directory separator character
%                 - cellstr
%                     { '\d{3}', '^PCASL.*', 'PCASL.*\.PAR$' } (i.e., each cell contains a regular expression to match a file or dir
%                 - structure array
%                       dirSpecs(1).pattern = '*';        % optional; only required to limit dir listing
%                       dirSpecs(1).regexp = '\d{3}';
%                       dirSpecs(2).pattern = 'PCASL_*';  % optional; only required to limit dir listing
%                       dirSpecs(2).regexp = '^PCASL.*';
%                       dirSpecs(3).pattern = '*.PAR';    % optional; only required to limit dir listing
%                       dirSpecs(3).regexp = '\.PAR$';
%     varargin  - Set of optional parameters provided in pairs of 'Name'-'Value':
%                   'Verbose'         - set to true or false to switch terminal feedback (DEFAULT FALSE)
%                   'IgnoreCase'      - set to true or false to ignore or respect case sensitivity during match (DEFAULT true on windows, false on linux)
%                   'StripRoot'       - do not include the root folder in the returned paths (DEFAULT FALSE)
%                   'FolderSeparator' - sets folder separator character(s) in dirspec string (DEFAULT '|')
%                   'Match'           - one of { 'Files', 'Directories', 'FilesAndDirectories' } (DEFAULT 'Files')
%                   'SortBy'          - vector containing token indices for sorting order; negative values sort descending
% OUTPUT:
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Recursively find files in the root directory according to the dirSpecs
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright © 2015-2019 ExploreASL

    % Setup argument parse for variable arguments
    p = inputParser;
    addOptional(p, 'Verbose', false, @islogical);
    addOptional(p, 'IgnoreCase', ispc(), @islogical);
    addOptional(p, 'StripRoot', false, @islogical);
    addOptional(p, 'FolderSeparator', '|', @ischar);
    possibleMatches = {'Files', 'Directories', 'FilesAndDirectories'};
    addOptional(p, 'Match', 'Files', @(x) any(validatestring(x,possibleMatches)));
    addOptional(p, 'SortBy', [], @isnumeric);

    % Parse parameters
    parse(p,varargin{:});
	
	if nargin < 2
		error('xasl_adm_FindByRegExp: Minimum number of required parameters is two.');
	end

    % Copy results to struct and convert some for easy use
    parms = p.Results;
    switch p.Results.Match
        case 'Files'
            parms.MatchFiles = true;
            parms.MatchDirectories = false;
        case 'Directories'
            parms.MatchFiles = false;
            parms.MatchDirectories = true;
        case 'FilesAndDirectories'
            parms.MatchFiles = true;
            parms.MatchDirectories = true;
    end

    % Check root folder
    if isempty(root)
        root = pwd;
    end
    if ~ischar(root)
        error('xasl_adm_FindByRegExp: First argument must be a string (root folder)');
    end
    if ~isdir(root)
        error('xasl_adm_FindByRegExp: Root folder doesn''t exist: %s', root);
    end
    
    % make sure the root ends with a separator and uses system native slashes
    root = fullfile(root,filesep); 
    if parms.StripRoot
        parms.nSkipChars = length(root);
    else
        parms.nSkipChars = 0;
    end
    
    %% check dirSpecs argument and allow several formats
    if isempty(dirSpecs)
        error('xasl_adm_FindByRegExp: dirSpecs argument must be a non-empty structure or string');
    end
    
    if ischar(dirSpecs)
        % accept strings like this:
        % 'root|regexp 1|regexp 1|regexp'
        dirSpecs = strsplit(dirSpecs, parms.FolderSeparator);
    end
    
    if iscellstr(dirSpecs)
        % accept cell string array's of regular expressions.
        newdirSpecs = struct('regexp',[]);
        for ii=1:length(dirSpecs)
            newdirSpecs(ii).regexp = dirSpecs{ii};
        end
        dirSpecs = newdirSpecs;
    end
    
    if ~isstruct(dirSpecs)
        error('xasl_adm_FindByRegExp: dirSpecs argument must be a string, cellstr or structure');
    end

    % Setup parameters for recursive calls
    currentlevel = 1;
    tree = {}; 
    mytokens = {};
    
    % and execute...
    [tree, mytokens] = FindByRegExpRec(root, dirSpecs, currentlevel, tree, parms, mytokens, {} );
    tree = tree'; % prefer rows for easy display
    
    % Sort output
    sortBy = parms.SortBy;
    if isempty(sortBy)
        sortBy = 1:size(mytokens,2);
    end
    [ mytokens, I] = sortrows(mytokens,sortBy);
    tree = tree(I);
    
    % Return tokens?
    if nargout>1
        optionalTokens = mytokens;
    end
end

% Internal recursive implementation
function [tree, rectokens] = FindByRegExpRec( root, dirSpecs, currentlevel, tree, parms, rectokens, parenttokens )

    L = dirSpecs(currentlevel);
	if parms.IgnoreCase
		casearg = 'ignorecase';
	else
		casearg = 'matchcase';
	end

	if isfield(L,'pattern')
		pattern = L.pattern;
	else
		pattern = '*';
	end

	entries = dir(fullfile(root, pattern));
	for iEntry=1:length(entries)
		E = entries(iEntry);
		if strcmp(E.name,'.') || strcmp(E.name,'..')
			continue; % ignore '.' and '..' (but accept hidden items if they are accepted by regexp)
		end
		
		bTryMatch = false;
		bRecurse = false;
		% scan folders if this is not the last level, otherwise handle filenames
		if currentlevel==length(dirSpecs)
			if (~E.isdir && parms.MatchFiles) || (E.isdir && parms.MatchDirectories)
				bTryMatch = true;
			end
		else
			if E.isdir
				bTryMatch = true;
				bRecurse = true;
			end
		end
		if bTryMatch
			[ match, tokens ] = regexp(E.name, L.regexp, 'match', 'tokens', 'once', casearg);
			
			% tricky: keep a copy of the parent tokens to be able to store the tokens for all matches
			if ~isempty(tokens)
				alltokens = [ parenttokens tokens ];
			else
				alltokens = parenttokens;
			end
			
			if ~isempty(match)
				% are we at the lowest matching directory level?
				if currentlevel==length(dirSpecs)
					% yes => add this match to the collection
					matchedpath = fullfile(root,E.name);
					if parms.nSkipChars>0
						matchedpath = matchedpath(parms.nSkipChars+1:end);
					end
					tree{end+1} = matchedpath; %#ok<AGROW>
					% tricky: copy the token off all directory levels (TODO: do we realy need a for to create a {N,LEVELS} array?)
					rectokens(length(tree),1:length(alltokens)) = alltokens;
					
					% be verbose?
					if parms.Verbose
						fprintf('match: %s\n',matchedpath);
					end
				end
				% enter sublevel if lowest level is not reached yet
				if bRecurse
					subroot = fullfile(root,E.name);
					[tree, rectokens] = FindByRegExpRec(subroot, dirSpecs, currentlevel+1, tree, parms, rectokens, alltokens);
				end
			end
		end
	end
end
