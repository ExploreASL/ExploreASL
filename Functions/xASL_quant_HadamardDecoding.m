function [imASLReordered] = xASL_quant_HadamardDecoding(imPath, xQ)
%xASL_quant_HadamardDecoding Hadamard-4 & Hadamard-8 Decoding
%
% FORMAT:       [imASLReordered] = xASL_quant_HadamardDecoding(imPath, xQ)
%
% INPUT:        imPath - Path to the encoded ASL4D image we want to decode (STRING, REQUIRED)
%
%               xQ     - x.Q field with Hadamard input parameters containing the following subfields
%                          - TimeEncodedMatrixType (REQUIRED)
%                             - Hadamard
%                             - Walsh
%                          - TimeEncodedMatrixSize (REQUIRED)
%                             - '4' for Hadamard-4
%                             - '8' for Hadamard-8
%                          - TimeEncodedMatrix (OPTIONAL)
%                             - Matrix given by the user
%                          - NumberEchoTimes  - Number of different echos (REQUIRED)
% OUTPUT:       imASLReordered - Decoded ASL volumes
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Hadamard-4 & Hadamard-8 Decoding.
%
% 0. Admin: Check inputs, load data, specify the decoding matrix
% 1. Reorder multi-TE data
% 2. Decode the Hadamard data
% 3. Reorder multi-TE back to the initial order of PLD/TE
% 4. Normalization of the decoded data
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      n/a
% __________________________________
% Copyright (c) 2015-2021 ExploreASL

%% 0. Admin
% Check if all inputs are present

if nargin<1 || isempty(imPath)
    error('imPath input is empty');
end

if nargin<2 || isempty(xQ)
    error('xQ input is empty');
end

% Specify the TimeEncodedMatrix
if isfield(xQ,'TimeEncodedMatrix') && ~isempty(xQ.TimeEncodedMatrix)
	if (size(xQ.TimeEncodedMatrix,1)+1) ~= size(xQ.TimeEncodedMatrix,2)
		error('TimeEncodedMatrix must be square (N-1)xN ');
	end
	
    % TimeEncodedMatrix exists, we verify the size
	if ~isfield(xQ,'TimeEncodedMatrixSize') || isempty(xQ.TimeEncodedMatrixSize)
		xQ.TimeEncodedMatrixSize = size(xQ.TimeEncodedMatrix,2);
	else
		if size(xQ.TimeEncodedMatrix,2) ~= xQ.TimeEncodedMatrixSize
			error('Mismatch between TimeEncodedMatrix and TimeEncodedMatrixSize');
		end
	end
else
    % TimeEncodedMatrix not provided, we must create it
	if ~isfield(xQ,'TimeEncodedMatrixType') || isempty(xQ.TimeEncodedMatrixType)
		error('Neither TimeEncodedMatrix and TimeEncodedMatrixType provided');
	end
	
	if ~isfield(xQ,'TimeEncodedMatrixSize') || isempty(xQ.TimeEncodedMatrixSize)
		error('Neither TimeEncodedMatrix and TimeEncodedMatrixSize provided');
	end
	
	% #### For Walsh Decoding Matrix ####
	if strcmp(xQ.TimeEncodedMatrixType,'Walsh')
		
		if xQ.TimeEncodedMatrixSize == 4
			xQ.TimeEncodedMatrix =...
				[1 -1  1 -1;
				 1 -1 -1  1;
				 1  1 -1 -1];
		elseif xQ.TimeEncodedMatrixSize == 8
            xQ.TimeEncodedMatrix = [1 -1  1 -1  1 -1  1 -1;
                                 1 -1  1 -1 -1  1 -1  1;
                                 1 -1 -1  1 -1  1  1 -1;
                                 1 -1 -1  1  1 -1 -1  1;
                                 1  1 -1 -1  1  1 -1 -1;
                                 1  1 -1 -1 -1 -1  1  1;
                                 1  1  1  1 -1 -1 -1 -1];
		end
		
		% #### For Hadamard Decoding Matrix ####
	elseif strcmp(xQ.TimeEncodedMatrixType,'Hadamard')
		if xQ.TimeEncodedMatrixSize == 4
            xQ.TimeEncodedMatrix = [1 -1  1 -1;
                                 1  1 -1 -1;
                                 1 -1 -1  1];
		elseif xQ.TimeEncodedMatrixSize == 8
            xQ.TimeEncodedMatrix = [1 -1 -1  1 -1  1  1 -1;
                                 1  1 -1 -1 -1 -1  1  1;
                                 1 -1  1 -1 -1  1 -1  1;
                                 1  1  1  1 -1 -1 -1 -1;
                                 1 -1 -1  1  1 -1 -1  1;
                                 1  1 -1 -1  1  1 -1 -1;
                                 1 -1  1 -1  1 -1  1 -1];
		end    
	end
end

% Load time-series nifti
imEncoded = xASL_io_Nifti2Im(imPath); 

% Calculate specific data sizes
nDecodedTI = xQ.TimeEncodedMatrixSize-1;                   % Number of TIs is always MatrixSize -1
nDecodedVolumes = xQ.NumberEchoTimes * nDecodedTI;

nEncodedVolumes = size(imEncoded, 4);

nPLDs = nEncodedVolumes / xQ.NumberEchoTimes;
nRepetitions = nEncodedVolumes / (xQ.TimeEncodedMatrixSize * xQ.NumberEchoTimes); % Calculating no. of acquisition repeats

%% 1. Reorder multi-TE data
% At this point the data is organized like this (in terms of ASL4D.nii volumes):
% PLD1/TE1,PLD1/TE2,PLD1/TE3,PLD1/TE4...PLD2/TE1,PLD2/TE2,PLD2/TE3... (PLDs first, TEs after)
%
% And for decoding we want
% TE1/PLD1,TE1/PLD2,TE1/PLD3,TE1/PLD4...TE2/PLD1,TE2/PLD2,TE2/PLD3,TE2/PLD4 (TEs first, PLDs after)
if isfield(xQ,'NumberEchoTimes') && (xQ.NumberEchoTimes > 1)
	% Reorder data - first cycle TEs afterwards PLDs
	vectorOldOrder = zeros(nEncodedVolumes, 1);
	for iTE = 1:xQ.NumberEchoTimes
		vectorOldOrder((1:nPLDs)+(iTE-1)*nPLDs) = (iTE-1)+1:xQ.NumberEchoTimes:nEncodedVolumes;
	end

	% Reorder the data
	imEncoded = imEncoded(:,:,:,vectorOldOrder);
end

%% 2. Decode the Hadamard data
Decoded_ASL = zeros(size(imEncoded,1), size(imEncoded,2), size(imEncoded,3), nDecodedVolumes * nRepetitions);    
idx=0;
for Repetition = 1:nRepetitions
    for TE = 1:xQ.NumberEchoTimes
        for TI = 1:nDecodedTI
            
            indexPositive = find(xQ.TimeEncodedMatrix(TI,:)==1);
            indexNegative = find(xQ.TimeEncodedMatrix(TI,:)==-1);
            Decoded_ASL(:,:,:,((TE-1)*nDecodedTI+TI)+(Repetition-1)*nDecodedVolumes) = mean(imEncoded(:,:,:,(indexPositive+idx)),4) - mean(imEncoded(:,:,:,(indexNegative+idx)),4);
            
        end
        idx = idx+xQ.TimeEncodedMatrixSize;
    end
end


%% 3. Reorder multi-TE back to the initial order of PLD/TE
%
% For model fitting, we want the PLDs-first-TEs-second order (just like the
% beginning) so we need to reorder it again

imASLReordered=xASL_im_HadamardReorder(Decoded_ASL,xQ.NumberEchoTimes);


%% 4. Normalization of the decoded data
% NormalizationFactor = 1/(m_INumSets/2);
% where m_INumSets is the number of images (e.g. 8 for Hadamard 8x8)

imASLReordered=xASL_im_HadamardDecodingNormalize(xQ.TimeEncodedMatrixSize, imASLReordered);

    
end


function imASLReordered=xASL_im_HadamardReorder(Decoded_ASL, NumberEchoTimes)
imEncoded = Decoded_ASL;
nPLDs = int32(size(imEncoded,4)/NumberEchoTimes); % this now takes size of decoded PLDs

vectorNewOrder = zeros(size(imEncoded,4),1);
for iPLD = 1:nPLDs
    vectorNewOrder((1:NumberEchoTimes)+(iPLD-1)*NumberEchoTimes) = (iPLD-1)+1:nPLDs:size(imEncoded,4);
end

imASLReordered = zeros(size(imEncoded,1),size(imEncoded, 2),size(imEncoded, 3),size(imEncoded, 4));
imASLReordered(:,:,:,1:end) = imEncoded(:,:,:,vectorNewOrder);

end

function imASLReordered=xASL_im_HadamardDecodingNormalize(TimeEncodedMatrixSize, imASLReordered)

NormalizationFactor = 1/(TimeEncodedMatrixSize/2);
imASLReordered = imASLReordered * NormalizationFactor;

end
