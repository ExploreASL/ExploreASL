function bTimeEncodedFME = xASL_imp_CheckIfFME(jsonIn, jsonOut, bTimeEncoded)
%xASL_imp_CheckIfFME Check if a sequence, based on its JSON, is a FME (Fraunhofer Mevis) time encoded sequence
%
% FORMAT: bTimeEncodedFME = xASL_imp_CheckIfFME(jsonIn[, jsonOut, bTimeEncoded])
% 
% INPUT:
%  jsonIn          - json structure of the sequence (REQUIRED)
%  jsonOut         - json structure of the sequence, already pre-parsed (OPTIONAL, DEFAULT=[])
%  bTimeEncoded    - true if known that the sequence is Time-Encoded (OPTIONAL, DEFAULT=0)
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
% Copyright 2015-2023 ExploreASL

    if nargin < 1 || isempty(jsonIn)
		error('Require JSON structure on the input');
	end
	
	if nargin < 2 || isempty(jsonOut)
		jsonOut = [];
	end
	
	if nargin < 3 || isempty(bTimeEncoded)
		bTimeEncoded = false;
	end
	
	% Determine if we have the specific FME Hadamard sequence from Bremen
	
	% This can be recognized either by the specific Series description
	bTimeEncodedFME = false;
	if isfield(jsonIn, 'SeriesDescription')
		if (~isempty(regexp(char(jsonIn.SeriesDescription),'(Encoded_Images_Had)\d\d(_)\d\d(_TIs_)\d\d(_TEs)', 'once')))% || ... % ASL format
				%~isempty(regexp(char(jsonIn.SeriesDescription),'(ss_TE)\d\d(_TI)\d\d\d\d', 'once')))                            % M0 format
			bTimeEncodedFME = true;
			return;
		end
	end
		
	% Or by the sequence name
	if isfield(jsonIn, 'SequenceName') && ~isempty(regexp(char(jsonIn.SequenceName), 'fme_asl', 'once'))
		if (isfield(jsonOut, 'PostLabelingDelay') && length(unique(jsonOut.PostLabelingDelay))>1)
			bTimeEncodedFME = true;
			return
		end
	
		if isfield(jsonOut, 'LabelingDuration') && length(unique(jsonOut.LabelingDuration))>1
			bTimeEncodedFME = true;
			return
		end
		
		if isfield(jsonOut, 'EchoTime') && length(unique(jsonOut.EchoTime))>1
			bTimeEncodedFME = true;
			return
		end
		
		if isfield(jsonOut, 'TimeEncodedMatrixSize') && ~isempty(jsonOut.TimeEncodedMatrixSize)
			bTimeEncodedFME = true;
			return
		end
		
		if bTimeEncoded
			bTimeEncodedFME = true;
			return
		end
				
	end

	% Or protocol name
	if isfield(jsonIn,'ProtocolName') && ~isempty(regexpi(char(jsonIn.ProtocolName), 'fme_gammastar_Had(4|8)', 'once'))
		bTimeEncodedFME = true;
	end
end