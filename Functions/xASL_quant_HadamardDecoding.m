function [imASL_reorder] = xASL_im_HadamardDecoding(Encoded_ASL, xDecodingFields, NumberEchoTimes)
%xASL_quant_HadamardDecodingHadamard-4 & Hadamard-8 Decoding
%
% FORMAT:       [imASL_reorder] = xASL_quant_HadamardDecoding(Encoded_ASL, xDecodingFields, NumberEchoTimes)
%
% INPUT:        Encoded_ASL      - ASL4D image we want to decode (REQUIRED)
%
%               xDecodingFields  1) TimeEncodedMatrixType (REQUIRED)
%                                   - Hadamard
%                                   - Walsh
%                                2) TimeEncodedMatrixSize (REQUIRED)
%                                   - '4' for Hadamard-4
%                                   - '8' for Hadamard-8
%                                3) DecodingMatrix (OPTIONAL)
%                                   - Matrix given by the user
%               
%               NumberEchoTimes  - Field from x.TimeEncodedEchoTimes (REQUIRED)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Hadamard-4 & Hadamard-8 Decoding.
%
% 1. Step-1: Reorder data
% 2. Step-2: decode data
% 3. Step-3: reorder again for model fitting
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      n/a
% __________________________________
% Copyright 2015-2021 ExploreASL

% Checking if all inputs are present

if nargin<3 && isempty(Encoded_ASL)
    warning('Encoded_ASL input is empty');
end

if nargin<3 && isempty(xDecodingFields)
    warning('xDecodingFields input is empty');
end

if nargin<3 && isempty(NumberEchoTimes)
    warning('NumberEchoTimes input is empty');
end

% Decoding Fields

TimeEncodedMatrixType = xDecodingFields.TimeEncodedMatrixType;
TimeEncodedMatrixSize = xDecodingFields.TimeEncodedMatrixSize;

if isfield(xDecodingFields,'DecodingMatrix') && ~isempty(xDecodingFields.DecodingMatrix)
    DecodingMatrix_input = xDecodingFields.DecodingMatrix;
else
    DecodingMatrix_input = [];
end

% #### For Walsh Decoding Matrix ####
if strcmp(TimeEncodedMatrixType,'Walsh')
    
    if TimeEncodedMatrixSize == 4
        
        ASL_im = xASL_io_Nifti2Im(Encoded_ASL); % Load time-series nifti
        
        if ~isempty(DecodingMatrix_input) %if there's a Decoding Matrix from the data
            DecodingMatrix = DecodingMatrix_input;
        else
            DecodingMatrix = [1 -1  1 -1;
                              1 -1 -1  1;
                              1  1 -1 -1];
        end
        
    elseif TimeEncodedMatrixSize == 8
        
        ASL_im = xASL_io_Nifti2Im(Encoded_ASL); % Load time-series nifti
        
        if ~isempty(DecodingMatrix_input) %if there's a Decoding Matrix from the data
            DecodingMatrix = DecodingMatrix_input;
        else
            DecodingMatrix = [1 -1  1 -1  1 -1  1 -1;
                              1 -1  1 -1 -1  1 -1  1;
                              1 -1 -1  1 -1  1  1 -1;
                              1 -1 -1  1  1 -1 -1  1;
                              1  1 -1 -1  1  1 -1 -1;
                              1  1 -1 -1 -1 -1  1  1;
                              1  1  1  1 -1 -1 -1 -1];
        end
    end
    
    % #### For Hadamard Decoding Matrix ####
elseif strcmp(TimeEncodedMatrixType,'Hadamard')
    
    
    
    if TimeEncodedMatrixSize == 4
        
        ASL_im = xASL_io_Nifti2Im(Encoded_ASL); % Load time-series nifti
        
        if ~isempty(DecodingMatrix_input) %if there's a Decoding Matrix from the data
            DecodingMatrix = DecodingMatrix_input;
        else
            DecodingMatrix = [1 -1  1 -1;
                              1  1 -1 -1;
                              1 -1 -1  1];
        end
        
    elseif TimeEncodedMatrixSize == 8
        
        ASL_im = xASL_io_Nifti2Im(Encoded_ASL); % Load time-series nifti
        
        if ~isempty(DecodingMatrix_input) %if there's a Decoding Matrix from the data
            DecodingMatrix = DecodingMatrix_input;
        else
            DecodingMatrix = [1 -1 -1  1 -1  1  1 -1;
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

nDecodedTI = TimeEncodedMatrixSize-1;                   % Number of TIs is always MatrixSize -1
nDecodedVolume = NumberEchoTimes * nDecodedTI;
nEncodedVolume = size(ASL_im, 4);
EncodedDataSize = TimeEncodedMatrixSize * NumberEchoTimes;     % Expected data size
nRepetitions = nEncodedVolume / EncodedDataSize;           % Calculating no. of acquisition repeats

DecodedDataSize = nDecodedVolume * nRepetitions;
Decoded_ASL = zeros(size(ASL_im,1),size(ASL_im,2),size(ASL_im,3),DecodedDataSize);
numberPLDs = int32(size(ASL_im,4)/NumberEchoTimes);

% Reorder data - first cycle TEs afterwards PLDs
vectorOldOrder = zeros(size(ASL_im,4),1);
for iTE = 1:(double(NumberEchoTimes))
    vectorOldOrder((1:numberPLDs)+(iTE-1)*numberPLDs) = (iTE-1)+1:NumberEchoTimes:size(ASL_im,4);
end
ASL_im(:,:,:,1:end) = ASL_im(:,:,:,vectorOldOrder);

%% step-2: decode data
    
idx=0;
for Repetition = 1:nRepetitions
    for TE = 1:NumberEchoTimes
        for TI = 1:nDecodedTI
            
            indexPositive = find(DecodingMatrix(TI,:)==1);
            indexNegative = find(DecodingMatrix(TI,:)==-1);
            Decoded_ASL(:,:,:,((TE-1)*nDecodedTI+TI)+(Repetition-1)*nDecodedVolume) = mean(ASL_im(:,:,:,(indexPositive+idx)),4) - mean(ASL_im(:,:,:,(indexNegative+idx)),4);
            
        end
        idx = idx+TimeEncodedMatrixSize;
    end
end


%% step-3: reorder again for model fitting
%
% For model fitting, we want the PLDs-first-TEs-second order (just like the
% beginning) so we need to reorder it again

imASL_reorder=xASL_im_HadamardReorder(Decoded_ASL,NumberEchoTimes);


%% ste-4: signal normalization
% NormalizationFactor = 1/(m_INumSets/2);
% where m_INumSets is the number of images (e.g. 8 for Hadamard 8x8)

imASL_reorder=xASL_im_HadamardDecodingNormalize(TimeEncodedMatrixSize,imASL_reorder);

    
end


function imASL_reorder=xASL_im_HadamardReorder(Decoded_ASL,NumberEchoTimes)
ASL_im = Decoded_ASL;
numberPLDs = int32(size(ASL_im,4)/NumberEchoTimes); % this now takes size of decoded PLDs

vectorNewOrder = zeros(size(ASL_im,4),1);
for iPLD = 1:numberPLDs
    vectorNewOrder((1:NumberEchoTimes)+(iPLD-1)*NumberEchoTimes) = (iPLD-1)+1:numberPLDs:size(ASL_im,4);
end

imASL_reorder = zeros(size(ASL_im,1),size(ASL_im,2),size(ASL_im,3),size(ASL_im,4));
imASL_reorder(:,:,:,1:end) = ASL_im(:,:,:,vectorNewOrder);

end

function imASL_reorder=xASL_im_HadamardDecodingNormalize(TimeEncodedMatrixSize,imASL_reorder)

NormalizationFactor = 1/(TimeEncodedMatrixSize/2);
imASL_reorder = imASL_reorder * NormalizationFactor;

end
