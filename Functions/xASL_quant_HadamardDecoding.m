function [imASLReordered] = xASL_quant_HadamardDecoding(imASLEncoded, xQ)
%xASL_quant_HadamardDecoding Hadamard-4 & Hadamard-8 Decoding
%
% FORMAT:       [imASLReordered] = xASL_quant_HadamardDecoding(imASLEncoded, xQ)
%
% INPUT:        imASLEncoded - ASL4D image we want to decode (REQUIRED)
%
%               xQ           - x.Q field with Hadamard input parameters containing the following subfields
%                              - TimeEncodedMatrixType (REQUIRED)
%                                   - Hadamard
%                                   - Walsh
%                              - TimeEncodedMatrixSize (REQUIRED)
%                                   - '4' for Hadamard-4
%                                   - '8' for Hadamard-8
%                              - TimeEncodedMatrix (OPTIONAL)
%                                   - Matrix given by the user
%                              - NumberEchoTimes  - Number of different echos (REQUIRED)
% OUTPUT:       imASLReordered - Decoded ASL volumes
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Hadamard-4 & Hadamard-8 Decoding.
%
% 1. Step-1: Reorder data
% 2. Step-2: Decode data
% 3. Step-3: Reorder again for model fitting
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      n/a
% __________________________________
% Copyright (c) 2015-2021 ExploreASL

% Check if all inputs are present

if nargin<1 || isempty(imASLEncoded)
    warning('imASLEncoded input is empty');
end

if nargin<2 || isempty(xQ)
    warning('xQ input is empty');
end

% Decoding Fields

if isfield(xQ,'TimeEncodedMatrix') && ~isempty(xQ.TimeEncodedMatrix)
    TimeEncodedMatrix_input = xQ.TimeEncodedMatrix;
else
    TimeEncodedMatrix_input = [];
end

% #### For Walsh Decoding Matrix ####
if strcmp(xQ.TimeEncodedMatrixType,'Walsh')
    
    if xQ.TimeEncodedMatrixSize == 4
        
        ASL_im = xASL_io_Nifti2Im(imASLEncoded); % Load time-series nifti
        
        if ~isempty(TimeEncodedMatrix_input) %if there's a Decoding Matrix from the data
            TimeEncodedMatrix = TimeEncodedMatrix_input;
        else
            TimeEncodedMatrix = [1 -1  1 -1;
                              1 -1 -1  1;
                              1  1 -1 -1];
        end
        
    elseif xQ.TimeEncodedMatrixSize == 8
        
        ASL_im = xASL_io_Nifti2Im(imASLEncoded); % Load time-series nifti
        
        if ~isempty(TimeEncodedMatrix_input) %if there's a Decoding Matrix from the data
            TimeEncodedMatrix = TimeEncodedMatrix_input;
        else
            TimeEncodedMatrix = [1 -1  1 -1  1 -1  1 -1;
                              1 -1  1 -1 -1  1 -1  1;
                              1 -1 -1  1 -1  1  1 -1;
                              1 -1 -1  1  1 -1 -1  1;
                              1  1 -1 -1  1  1 -1 -1;
                              1  1 -1 -1 -1 -1  1  1;
                              1  1  1  1 -1 -1 -1 -1];
        end
    end
    
    % #### For Hadamard Decoding Matrix ####
elseif strcmp(xQ.TimeEncodedMatrixType,'Hadamard')
    
    
    
    if xQ.TimeEncodedMatrixSize == 4
        
        ASL_im = xASL_io_Nifti2Im(imASLEncoded); % Load time-series nifti
        
        if ~isempty(TimeEncodedMatrix_input) %if there's a Decoding Matrix from the data
            TimeEncodedMatrix = TimeEncodedMatrix_input;
        else
            TimeEncodedMatrix = [1 -1  1 -1;
                              1  1 -1 -1;
                              1 -1 -1  1];
        end
        
    elseif xQ.TimeEncodedMatrixSize == 8
        
        ASL_im = xASL_io_Nifti2Im(imASLEncoded); % Load time-series nifti
        
        if ~isempty(TimeEncodedMatrix_input) %if there's a Decoding Matrix from the data
            TimeEncodedMatrix = TimeEncodedMatrix_input;
        else
            TimeEncodedMatrix = [1 -1 -1  1 -1  1  1 -1;
                              1  1 -1 -1 -1 -1  1  1;
                              1 -1  1 -1 -1  1 -1  1;
                              1  1  1  1 -1 -1 -1 -1;
                              1 -1 -1  1  1 -1 -1  1;
                              1  1 -1 -1  1  1 -1 -1;
                              1 -1  1 -1  1 -1  1 -1];
        end
    end
    
end

%% step-1: Reorder data
% At this point the data is organized like this (in terms of ASL4D.nii volumes):
% PLD1/TE1,PLD1/TE2,PLD1/TE3,PLD1/TE4...PLD2/TE1,PLD2/TE2,PLD2/TE3... (PLDs first, TEs after)
%
% And for decoding we want
% TE1/PLD1,TE1/PLD2,TE1/PLD3,TE1.PL4...TE2/PLD1,TE2/PLD2,TE2/PLD3,TE2/PLD4 (TEs first, PLDs after)

nDecodedTI = xQ.TimeEncodedMatrixSize-1;                   % Number of TIs is always MatrixSize -1
nDecodedVolume = xQ.NumberEchoTimes * nDecodedTI;
nEncodedVolume = size(ASL_im, 4);
EncodedDataSize = xQ.TimeEncodedMatrixSize * xQ.NumberEchoTimes;     % Expected data size
nRepetitions = nEncodedVolume / EncodedDataSize;           % Calculating no. of acquisition repeats

DecodedDataSize = nDecodedVolume * nRepetitions;
Decoded_ASL = zeros(size(ASL_im,1),size(ASL_im,2),size(ASL_im,3),DecodedDataSize);
numberPLDs = int32(size(ASL_im,4)/xQ.NumberEchoTimes);

% Reorder data - first cycle TEs afterwards PLDs
vectorOldOrder = zeros(size(ASL_im,4),1);
for iTE = 1:(double(xQ.NumberEchoTimes))
    vectorOldOrder((1:numberPLDs)+(iTE-1)*numberPLDs) = (iTE-1)+1:xQ.NumberEchoTimes:size(ASL_im,4);
end
ASL_im(:,:,:,1:end) = ASL_im(:,:,:,vectorOldOrder);

%% step-2: decode data
    
idx=0;
for Repetition = 1:nRepetitions
    for TE = 1:xQ.NumberEchoTimes
        for TI = 1:nDecodedTI
            
            indexPositive = find(TimeEncodedMatrix(TI,:)==1);
            indexNegative = find(TimeEncodedMatrix(TI,:)==-1);
            Decoded_ASL(:,:,:,((TE-1)*nDecodedTI+TI)+(Repetition-1)*nDecodedVolume) = mean(ASL_im(:,:,:,(indexPositive+idx)),4) - mean(ASL_im(:,:,:,(indexNegative+idx)),4);
            
        end
        idx = idx+xQ.TimeEncodedMatrixSize;
    end
end


%% step-3: reorder again for model fitting
%
% For model fitting, we want the PLDs-first-TEs-second order (just like the
% beginning) so we need to reorder it again

imASLReordered=xASL_im_HadamardReorder(Decoded_ASL,xQ.NumberEchoTimes);


%% ste-4: signal normalization
% NormalizationFactor = 1/(m_INumSets/2);
% where m_INumSets is the number of images (e.g. 8 for Hadamard 8x8)

imASLReordered=xASL_im_HadamardDecodingNormalize(xQ.TimeEncodedMatrixSize, imASLReordered);

    
end


function imASLReordered=xASL_im_HadamardReorder(Decoded_ASL, NumberEchoTimes)
ASL_im = Decoded_ASL;
numberPLDs = int32(size(ASL_im,4)/NumberEchoTimes); % this now takes size of decoded PLDs

vectorNewOrder = zeros(size(ASL_im,4),1);
for iPLD = 1:numberPLDs
    vectorNewOrder((1:NumberEchoTimes)+(iPLD-1)*NumberEchoTimes) = (iPLD-1)+1:numberPLDs:size(ASL_im,4);
end

imASLReordered = zeros(size(ASL_im,1),size(ASL_im, 2),size(ASL_im, 3),size(ASL_im, 4));
imASLReordered(:,:,:,1:end) = ASL_im(:,:,:,vectorNewOrder);

end

function imASLReordered=xASL_im_HadamardDecodingNormalize(TimeEncodedMatrixSize, imASLReordered)

NormalizationFactor = 1/(TimeEncodedMatrixSize/2);
imASLReordered = imASLReordered * NormalizationFactor;

end
