function strOut = xASL_adm_CorrectName(strIn, bOption, strExclude)
% Go through a string and remove all non-word characters
%
% FORMAT: strOut = xASL_adm_CorrectName(strIn[, bOption, strExclude])
%
% INPUT:
%   strIn  	   - input string
%   bOption    - how to deal with the non-word characters (DEFAULT 1)
%                1 replaces non-word characters with underscore, removes double underscores
%                2 replaces non-word characters (including underscore) by empty space, concatenating all words
%   strExclude - list of characters that are exempt from the replacement (DEFAULT '')
%
% OUTPUT:
%   strOut     - corrected output string
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Finds and replaces all non-word characters either by empty space or by an underscore. 
%              Optionally leaves in few selected special characters. Note that if '_' is excluded from 
%              replacement, but option 2 is on, then underscores are replaced anyway.
%
% EXAMPLE:     xASL_adm_CorrectName('t_es._-ting^-&&@function')         returns 't_es_ting_function'
%              xASL_adm_CorrectName('t_es._-ting^-&&@function',1)       returns 't_es_ting_function'
%              xASL_adm_CorrectName('t_es._-ting^-&&@function',2)       returns 'testingfunction'
%              xASL_adm_CorrectName('t_es._-ting^-&&@function',1,'&-')  returns 't_es_-ting_-&&_function'
%              xASL_adm_CorrectName('t_es._-ting^-&&@function',2,'&-')  returns 'tes-ting-&&function'
%              xASL_adm_CorrectName('t_es._-ting^-&&@function',2,'&-_') returns 't_es_-ting-&&function'
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright ?? 2015-2020 ExploreASL
%
% 2019-04-05 JP

    % Admin
	% By default use Option 1
	if nargin < 2 || isempty(bOption)
		bOption = 1;
	end
	
	% By default, exclude no strings
	if nargin < 3
		strExclude = [];
	end
	
	% Don't take more than three input arguments
	if nargin > 3
		error('xASL_adm_CorrectName: Too many input parameters.');
    end
	
    if ~ischar(strIn)
        warning('Input string wasnt a string, check input format');
        strOut = strIn;        
        return;
    end
    
	% Making own list is deprecated - we use the regexp function now
	IllegalList = {' ' '*' '.' '"' '/' '\' '[' ']' '{' '}' ':' ';' '+' '=' ',' '<' '>' '?' '~' '@' '#' '$' '%' '^' '&' '*' '(' ')' '-' '|'};
	%for iI=1:length(IllegalList)
	% IndicesN            = find(strIn==IllegalList{iI});
	%	strOut(IndicesN)     = '_';
	%end
	
    % Here we define a character to replace the illegal characters with
    % By default this is an underscore, but if bOption==2 and we want to
    % keep the underscore (i.e. inside strExclude) we need something else
    
    ReplaceString = '_'; % default
    if bOption==2 && ~isempty(strExclude) && ~isempty(regexp(strExclude,'_'))
        IndicesString = find(cellfun(@(y) isempty(regexp(strExclude,y)), IllegalList(2:end)));
        if ~isempty(IndicesString)
            ReplaceString = IllegalList{IndicesString(1)+1};
        end % else keep default
    end
    strOut(1:length(strIn)) = ReplaceString;
	
	% Find the wrong characters
	strReg = '\w';
	if ~isempty(strExclude)
		strReg(4:2:(4+(length(strExclude)-1)*2)) = strExclude;
		strReg(3:2:(3+(length(strExclude)-1)*2)) = '|';
	end
	
	indPass = regexp(strIn,strReg);
	
	strOut(indPass) = strIn(indPass);
	
	if bOption==2
		% Remove all replacements
		strOut = strOut(strOut~=ReplaceString);
	else
		% Correct for double underscores
		[Ind1, Ind2] = regexp(strOut,'_*','start','end');
		Indices = true(1,length(strOut));
		if length(Ind1)==length(Ind2) % QA
			for iI=1:length(Ind1)
				Indices(Ind1(iI):Ind2(iI)-1) = false;
			end
		end
		strOut = strOut(Indices);
	end
end

