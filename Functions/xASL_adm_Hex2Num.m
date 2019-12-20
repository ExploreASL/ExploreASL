function outNum = xASL_adm_Hex2Num(inStr, type, endian)
% Converts string in hexadecimal to a decimal number. 
%
% FORMAT: outNum = xASL_adm_hex2num(inStr)
%
% INPUT:
%   inStr     - string containing the number in hexadecimal format
%   type      - 'FLOAT' (16,32,64,128-bit) floating point IEEE 754:1985
%               'DECIMAL' (normally written - ANSI X3.9)
%               'SINT'
%               'UINT'
%   endian    - 0 little endian (DEFAULT) 
%               1 big endian (DEFAULT for DECIMAL)
% OUTPUT:
%   outNum    - converted to decimal float/double
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Takes a hexadecimal string and converts it to number. Works also when the string contains escape characters, and for single-floats and 
%              for a little and big endian
%              If containing 8 and less characters than treat as float, if more than as double
%
% EXAMPLE: outNum = xASL_adm_hex2num('b1bd5039')
%          outNum = xASL_adm_hex2num('b1/bd/50/39')
%          outNum = xASL_adm_hex2num('b1bd50392040')
%          outNum = xASL_adm_hex2num('3130303132322e31',2)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright © 2015-2019 ExploreASL
%
% 2019-03-29 JP

% Not enough input parameters
if (nargin<1) || isempty(inStr)
	outNum = NaN;
	return;
end

if (nargin<2) || isempty(type)
	type = 'float';
end

if (nargin<3) || isempty(endian)
	if strcmpi(type,'decimal')
		endian = 1;
	else
		endian = 0;
	end
end

% Takes only char on input
if ~ischar(inStr)
	outNum = NaN;
	return;
end

% Take out the backslashes
outNum = xASL_adm_CorrectName(inStr,2);

% Skip padding with zeros for now...
% 
% 	% Pad with zeros to 8 characters
% 	if length(outNum) < 8
% 		if endian
% 			outNum((end+1):8) = '0';
% 		else
% 			tmp = outNum;
% 			outNum = repmat('0',[1 8]);
% 			outNum((8-length(tmp)+1):8) = tmp;
% 		end
% 	end
% 		if length(outNum) < 32
% 			if endian
% 				outNum((end+1):32) = '0';
% 			else
% 				tmp = outNum;
% 				outNum = repmat('0',[1 32]);
% 				outNum((32-length(outNum)+1):32) = tmp;
% 			end
% 		end

% Needs to have even length
if mod(length(outNum),2)
	error('xASL_adm_Hex2Num: Even number of characters is required for a HEX number.');
end
	
% Flip if little endian
if ~endian
	tmp = outNum;
	for ii=1:(length(outNum)/2)
		outNum((1:2)+(2*(ii-1))) = tmp((length(outNum)+[-1,0]) - 2*(ii-1));
	end
end

switch (lower(type))
	case 'decimal'
		% Delivered as a decimal string in hex
	
		% First convert hex to decimal
		outNum = hex2dec(reshape(outNum,2,[])');
		
		% Then decimal to char
		outNum = char(outNum');
		
		% then convert char to nubmer
		outNum = str2double(outNum);
		
	case 'float'
		switch(length(outNum))
			case 4
				% 16 bit
				% Otherwise use a double conversion
				outBin = dec2bin(hex2dec(outNum),16);
				
				% First bit is sign
				if outBin(1) == '1'
					outNum = -1;
				else
					outNum = 1;
				end
				
				% Bits 2:9 gives the exponent
				outNum = outNum * (2^(bin2dec(outBin(2:6))-15));
				
				% Bits 10:32 give the mantissa
				% Converts to binary
				mantissa = (outBin(7:16) == '1');
				mantissa = 1+sum(mantissa.* (2.^(-1:-1:-10)));
				outNum = outNum * mantissa;
			case 8
				% It is a single precision - we need to do the conversion ourselves
				outBin = dec2bin(hex2dec(outNum),32);
				
				% First bit is sign
				if outBin(1) == '1'
					outNum = -1;
				else
					outNum = 1;
				end
				
				% Bits 2:9 gives the exponent
				outNum = outNum * (2^(bin2dec(outBin(2:9))-127));
				
				% Bits 10:32 give the mantissa
				% Converts to binary
				mantissa = (outBin(10:32) == '1');
				mantissa = 1+sum(mantissa.* (2.^(-1:-1:-23)));
				outNum = outNum * mantissa;
			case 16
				% Otherwise use a double conversion
				outBin = dec2bin(hex2dec(outNum),64);
				
				% First bit is sign
				if outBin(1) == '1'
					outNum = -1;
				else
					outNum = 1;
				end
				
				% Bits 2:9 gives the exponent
				outNum = outNum * (2^(bin2dec(outBin(2:12))-1023));
				
				% Bits 10:32 give the mantissa
				% Converts to binary
				mantissa = (outBin(13:64) == '1');
				mantissa = 1+sum(mantissa.* (2.^(-1:-1:-52)));
				outNum = outNum * mantissa;
			case 32
				% 128 bit
				% Otherwise use a double conversion
				outBin = dec2bin(hex2dec(outNum),64);
				
				% First bit is sign
				if outBin(1) == '1'
					outNum = -1;
				else
					outNum = 1;
				end
				
				% 15 Bits gives the exponent
				outNum = outNum * (2^(bin2dec(outBin(2:16))-16383));
				
				% Bits 17:128 give the mantissa
				% Converts to binary
				mantissa = (outBin(17:128) == '1');
				mantissa = 1+sum(mantissa.* (2.^(-1:-1:-112)));
				outNum = outNum * mantissa;
			otherwise
				error('xASL_adm_Hex2Num: Floating point conversion done only for 16,32,64,128-bit.');
		end
	case 'sint'
		outBin = dec2bin(hex2dec(outNum),64);
		% First bit is sign
		if outBin(1) == '1'
			outNum = -1;
		else
			outNum = 1;
		end
		outNum = outNum*sum((outBin(2:end)=='1').*(2.^((length(outBin)-2):-1:0)));
	case 'uint'
		outBin = dec2bin(hex2dec(outNum),64);
		outNum = sum((outBin(1:end)=='1').*(2.^((length(outBin)-1):-1:0)));
	otherwise
		error(['xASL_adm_Hex2Num: Unknown type ' type]);
end

