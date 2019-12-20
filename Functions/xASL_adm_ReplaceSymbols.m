function strOut = xASL_adm_ReplaceSymbols(strIn, symbolTable, bracketLeft, bracketRight)
% Helper function to replace <symbols> in strings using the values in the symbol table.
%
% FORMAT: 
%     strOut = xASL_adm_ReplaceSymbols(strIn, symbolTable[, bracketLeft, bracketRight])
% INPUT:
%     strIn        - input string (REQUIRED)
%     symbolTable  - table of symbols (REQUIRED)
%     bracketLeft  - string describing the left side bracket (OPTIONAL, DEFAULT = '<')
%     bracketright - string describing the right side bracket (OPTIONAL, DEFAULT = '<')
% OUTPUT:
%     strOut - replaced string
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: It takes the STRIN on input, then looks for symbols between BRACKETLEFT and BRACKETRIGHT and replaces these symbols in
%              in the string by the values provided in the SYMBOLTABLE as SYMBOLTABLE.SYMBOL, SYMBOLTABLE.D.SYMBOL, or SYMBOLTABLE.P.SYMBOL
%
% EXAMPLE: symbolTable.symA   = 'AAA';
%          symbolTable.symB   = 'BBB'; 
%          symbolTable.D.symC = 'CCC'; 
%          strOut = xASL_adm_ReplaceSymbols('Test<symA>-<symB>-<symC>, symbolTable,'<','>');
%          strOut is 'TestAAA-BBB-CCC'.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright © 2015-2019 ExploreASL

    % Admin
	if ~ischar(strIn)
		error('xASL_adm_ReplaceSymbols only works on text strings, not on %s', class(strIn));
	end
	
	if nargin<2 || isempty(symbolTable)
		error('xASL_adm_ReplaceSymbols needs at least two input parameters.');
	end

	if nargin<3 || isempty(bracketLeft)
		bracketLeft='<';
	end
	
	if nargin<4 || isempty(bracketRight)
		bracketRight='>';
	end
	
    bracketLeftLen = length(bracketLeft);
    bracketLeft = regexptranslate('escape',bracketLeft);
    
    bracketRightLen = length(bracketRight);
    bracketRight = regexptranslate('escape',bracketRight);
    
    max_depth = 10;
	strOut = strIn;
    for iRecurse=1:max_depth
        symbols = regexp(strOut, [bracketLeft '[a-zA-Z]\w*' bracketRight], 'match');
        if isempty(symbols)
            break; % no symbols anymore
        elseif iRecurse==max_depth
            error('xASL_adm_ReplaceSymbols: Symbol definition is probably cyclic: %s\nAborted expansion after %u iterations!\n',strOut,max_depth)
        end
        for iSymbol=1:length(symbols)
            bracket = symbols{iSymbol}; % ie '<SYMBOLNAME>'
            symbolName = bracket(bracketLeftLen+1:end-bracketRightLen); % ie 'SYMBOLNAME'
            if      isfield(symbolTable,symbolName)
                    symbolValue = symbolTable.(symbolName);
            elseif  isfield(symbolTable.D,symbolName) % directories
                    symbolValue = symbolTable.D.(symbolName);
            elseif  isfield(symbolTable.P,symbolName) % paths
                    symbolValue = symbolTable.P.(symbolName);
            else
                    error('xASL_adm_ReplaceSymbols: Symbol not defined: %s',symbolName);
            end
                
            if isnumeric(symbolValue)
                symbolValue = num2str(symbolValue);
            elseif ~ischar(symbolValue)
                error('xASL_adm_ReplaceSymbols: Cannot insert symbols that are not text strings: %s is a %s', symbolName, class(symbolValue))
            end
            strOut = strrep(strOut,bracket,symbolValue);
        end
    end
end
