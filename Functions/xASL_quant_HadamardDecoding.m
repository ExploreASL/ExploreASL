function xASL_quant_HadamardDecoding(Encoded_ASL, xDecodingFields, NumberEchoTimes)
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
    
    %.... Hadamard (natural) matrices....
    
end

%% step-1: Reorder data

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
ASL_im = Decoded_ASL;
numberPLDs = int32(size(ASL_im,4)/numberTE); % this now takes size of decoded PLDs

vectorNewOrder = zeros(size(ASL_im,4),1);
for iPLD = 1:(double(numberPLDs))
    vectorNewOrder((1:numberTE)+(iPLD-1)*numberTE) = (iPLD-1)+1:numberPLDs:size(ASL_im,4);
end

imASL_reorder = zeros(size(ASL_im,1),size(ASL_im,2),size(ASL_im,3),size(ASL_im,4));
imASL_reorder(:,:,:,1:end) = ASL_im(:,:,:,vectorNewOrder);

% Save 4D nifti file
xASL_io_SaveNifti(Encoded_ASL, Encoded_ASL, imASL_reorder)

    
    % For visualization/testing -> remove this afterwards
%     figure 
%     image(Decoded_ASL(:,:,12,5));
    
end