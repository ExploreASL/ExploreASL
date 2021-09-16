function [imASL_reorder] = xASL_quant_HadamardDecoding(Encoded_ASL, xDecodingFields, NumberEchoTimes)
%xASL_im_HadamardDecoding Hadamard-4 & Hadamard-8 Decoding
%
% FORMAT:       [Decoded_ASL] = xASL_im_HadamardDecoding(HadamardType, sec)
%
% INPUT:        HadamardType  - 4 for Hadamard4
%                             - 8 for Hadamard8
%                             - 12 for Hadamard12
%               sec           - output from xASL_im_HadamardCaseByCase
%
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Hadamard-4 & Hadamard-8 Decoding.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      n/a
% __________________________________
% Copyright 2015-2021 ExploreASL

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

numberTE = NumberEchoTimes;
nDecodedTI = TimeEncodedMatrixSize-1;                   % Number of TIs is always MatrixSize -1
nDecodedVol = numberTE * nDecodedTI;
nEncodedVol = size(ASL_im, 4);
EncodedDataSize = TimeEncodedMatrixSize * numberTE;     % Expected data size
nRepetitions = nEncodedVol / EncodedDataSize;           % Calculating no. of acquisition repeats

DecodedDataSize = nDecodedVol * nRepetitions;
Decoded_ASL = zeros(size(ASL_im,1),size(ASL_im,2),size(ASL_im,3),DecodedDataSize);
numberPLDs = int32(size(ASL_im,4)/numberTE);

% Reorder data - first cycle TEs afterwards PLDs
vectorOldOrder = zeros(size(ASL_im,4),1);
for iTE = 1:(double(numberTE))
    vectorOldOrder((1:numberPLDs)+(iTE-1)*numberPLDs) = (iTE-1)+1:numberTE:size(ASL_im,4);
end
ASL_im(:,:,:,1:end) = ASL_im(:,:,:,vectorOldOrder);

%% step-2: decode data
    
idx=0;
for Repetition = 1:nRepetitions
    for TE = 1:numberTE
        for TI = 1:nDecodedTI
            
            indexPositive = find(DecodingMatrix(TI,:)==1);
            indexNegative = find(DecodingMatrix(TI,:)==-1);
            Decoded_ASL(:,:,:,((TE-1)*nDecodedTI+TI)+(Repetition-1)*nDecodedVol) = mean(ASL_im(:,:,:,(indexPositive+idx)),4) - mean(ASL_im(:,:,:,(indexNegative+idx)),4);
            
        end
        idx=idx+TimeEncodedMatrixSize;
    end
end


%% step-3: reorder again for model fitting
%
% For model fitting, we want the PLDs-first-TEs-second order (just like the
% beginning) so we need to reorder it again

ASL_im = Decoded_ASL;
numberPLDs = int32(size(ASL_im,4)/numberTE); % this now takes size of decoded PLDs

vectorNewOrder = zeros(size(ASL_im,4),1);
for iPLD = 1:(double(numberPLDs))
    vectorNewOrder((1:numberTE)+(iPLD-1)*numberTE) = (iPLD-1)+1:numberPLDs:size(ASL_im,4);
end

imASL_reorder = zeros(size(ASL_im,1),size(ASL_im,2),size(ASL_im,3),size(ASL_im,4));
imASL_reorder(:,:,:,1:end) = ASL_im(:,:,:,vectorNewOrder);

%% ste-4: signal normalization
% NormalizationFactor = 1/(m_INumSets/2);
% where m_INumSets is the number of images (e.g. 8 for Hadamard 8x8)

NormalizationFactor = 1/(TimeEncodedMatrixSize/2);
imASL_reorder = imASL_reorder * NormalizationFactor;

    
end