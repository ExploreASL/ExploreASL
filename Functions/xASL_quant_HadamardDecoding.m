function [imDecoded] = xASL_quant_HadamardDecoding(imPath, xQ)
%xASL_quant_HadamardDecoding Hadamard-4 & Hadamard-8 Decoding
%
% FORMAT:       [imDecoded] = xASL_quant_HadamardDecoding(imPath, xQ)
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
% OUTPUT:       imDecoded - Decoded ASL volumes
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Hadamard-4 & Hadamard-8 Decoding.
%
% 0. Admin: Check inputs, load data
% 1. Specify the decoding matrix
% 2. Reorder multi-TE data
% 3. Decode the Hadamard data
% 4. Reorder multi-TE back to the initial order of PLD/TE
% 5. Normalization of the decoded data
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

%% 1. Specify the Encoding matrix

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
    
    if strcmp(xQ.Vendor,'Philips')
        xQ.TimeEncodedMatrix = xQ.TimeEncodedMatrix * -1;
    end
    
end

% Load time-series nifti
imEncoded = xASL_io_Nifti2Im(imPath); 

% Calculate specific data sizes
nDecodedTI = xQ.TimeEncodedMatrixSize-1;                   % Number of TIs is always MatrixSize -1
nDecodedVolumes = xQ.NumberEchoTimes * nDecodedTI;

nEncodedVolumes = size(imEncoded, 4);

nRepetitions = nEncodedVolumes / (xQ.TimeEncodedMatrixSize * xQ.NumberEchoTimes); % Calculating no. of acquisition repeats

%% 2. Reorder multi-TE data
% At this point the data is organized like this (in terms of ASL4D.nii volumes):
% PLD1/TE1,PLD1/TE2,PLD1/TE3,PLD1/TE4...PLD2/TE1,PLD2/TE2,PLD2/TE3... (TEs in first dimension, PLDs after)
%
% And for decoding we want
% TE1/PLD1,TE1/PLD2,TE1/PLD3,TE1/PLD4...TE2/PLD1,TE2/PLD2,TE2/PLD3,TE2/PLD4 (PLDs in the first dimension, TEs after)
if isfield(xQ,'NumberEchoTimes') && (xQ.NumberEchoTimes > 1)
	% This is the original order
	vectorOldOrder = 1:nEncodedVolumes;
	
	% Shape to a matrix with TEs first and all the rest later
	vectorOldOrder = reshape(vectorOldOrder, xQ.NumberEchoTimes, nEncodedVolumes/xQ.NumberEchoTimes);
	
	% Flip the two dimensions and make a row vector again
	vectorOldOrder = reshape(vectorOldOrder', 1, nEncodedVolumes);

	% Reorder the data
	imEncoded = imEncoded(:,:,:,vectorOldOrder);
end

%% 3. Decode the Hadamard data
imDecoded = zeros(size(imEncoded,1), size(imEncoded,2), size(imEncoded,3), nDecodedVolumes * nRepetitions);    
idx=0;
for iRepetition = 1:nRepetitions
    for iTE = 1:xQ.NumberEchoTimes
        for iTI = 1:nDecodedTI
            
            indexPositive = find(xQ.TimeEncodedMatrix(iTI,:)==1);
            indexNegative = find(xQ.TimeEncodedMatrix(iTI,:)==-1);
            imDecoded(:,:,:,((iTE-1)*nDecodedTI+iTI)+(iRepetition-1)*nDecodedVolumes) = mean(imEncoded(:,:,:,(indexPositive+idx)),4) - mean(imEncoded(:,:,:,(indexNegative+idx)),4);
            
        end
        idx = idx+xQ.TimeEncodedMatrixSize;
    end
end


%% 4. Reorder multi-TE back to the initial order of PLD/TE
%
% For model fitting, we want the PLDs-first-TEs-second order (just like the
% beginning) so we need to reorder it again
if isfield(xQ,'NumberEchoTimes') && (xQ.NumberEchoTimes > 1)
	% This is the original order
	vectorOldOrder = 1:size(imDecoded, 4);
	
	% Shape to a matrix with TEs second and all the rest first
	vectorOldOrder = reshape(vectorOldOrder, size(imDecoded, 4)/xQ.NumberEchoTimes, xQ.NumberEchoTimes);
	
	% Flip the two dimensions and make a row vector again
	vectorOldOrder = reshape(vectorOldOrder', 1, size(imDecoded, 4));

	% Reorder the data
	imDecoded = imDecoded(:,:,:,vectorOldOrder);
end

%% 5. Normalization of the decoded data
% NormalizationFactor = 1/(m_INumSets/2);
% where m_INumSets is the number of images (e.g. 8 for Hadamard 8x8)

imDecoded = imDecoded / (xQ.TimeEncodedMatrixSize/2);
    
end
