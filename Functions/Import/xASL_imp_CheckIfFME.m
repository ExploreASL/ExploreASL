function bTimeEncodedFME = xASL_imp_CheckIfFME(jsonIn)
%xASL_imp_CheckIfFME Check if a sequence, based on its JSON, is a FME (Fraunhofer Mevis) time encoded sequence
%
% FORMAT: bTimeEncodedFME = xASL_imp_CheckIfFME(jsonIn)
% 
% INPUT:
%  jsonIn          - json structure of the sequence (REQUIRED)
%
% OUTPUT:
%   bTimeEncodedFME    - Time encoded FME sequence or not (BOOLEAN)
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Check if the sequence, based on its JSON structure, is a FME (Fraunhofer Mevis) time encoded sequence.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     n/a
%
% __________________________________
% Copyright 2015-2022 ExploreASL

    if nargin < 1 || isempty(jsonIn)
		error('Require JSON structure on the input');
	end
	
	% Determine if we have the specific FME Hadamard sequence from Bremen
	
	% This can be recognized either by the specific Series description
	if isfield(jsonIn, 'SeriesDescription') &&...
			(~isempty(regexp(char(jsonIn.SeriesDescription),'(Encoded_Images_Had)\d\d(_)\d\d(_TIs_)\d\d(_TEs)', 'once')) || ... % ASL format
			~isempty(regexp(char(jsonIn.SeriesDescription),'(ss_TE)\d\d(_TI)\d\d\d\d', 'once')))                            % M0 format
		bTimeEncodedFME = true;
	% Or by the sequence name
	elseif isfield(jsonIn, 'SequenceName') && ~isempty(regexp(char(jsonIn.SequenceName), 'fme_asl', 'once'))
		bTimeEncodedFME = true;
	else
		bTimeEncodedFME = false;
	end
end