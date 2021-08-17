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

% #### For Hadamard Decoding Matrix ####

if strcmp(TimeEncodedMatrixType,'Hadamard')
    
    if TimeEncodedMatrixSize == 4
        
        ASL_im = xASL_io_Nifti2Im(Encoded_ASL); % Load time-series nifti
        
        if ~isempty(DecodingMatrix_input) %if there's a Decoding Matrix from the data
            DecodingMatrix = DecodingMatrix_input;
        else
            DecodingMatrix = [1 -1 1 -1;
                1 -1 -1 1;
                1 1  -1 -1];
        end
        
    elseif TimeEncodedMatrixSize == 8
        
        ASL_im = xASL_io_Nifti2Im(Encoded_ASL); % Load time-series nifti
        
        if ~isempty(DecodingMatrix_input) %if there's a Decoding Matrix from the data
            DecodingMatrix = DecodingMatrix_input;
        else
            DecodingMatrix = [1 -1 1 -1 1 -1 1 -1;
                1 -1 1 -1 -1 1 -1 1;
                1 -1 -1 1 -1 1 1 -1;
                1 -1 -1 1 1 -1 -1 1;
                1 1 -1 -1 1 1 -1 -1;
                1 1 -1 -1 -1 -1 1 1;
                1 1 1 1 -1 -1 -1 -1];
        end
        
    end
    
% #### For Walsh Decoding Matrix ####
elseif strcmp(TimeEncodedMatrixType,'Walsh')
    
    
    
    %.... Walsh matrices....
    
end

    numberTE = NumberEchoTimes;                           
    nDecodedTI = TimeEncodedMatrixSize-1;                   % Number of TIs is always MatrixSize -1
    nDecodedVol = numberTE * nDecodedTI;
    nEncodedVol = size(ASL_im, 4);                          
    EncodedDataSize = TimeEncodedMatrixSize * numberTE;     % Expected data size
    nRepetitions = nEncodedVol / EncodedDataSize;           % Calculating no. of acquisition repeats 
    
    DecodedDataSize = nDecodedVol * nRepetitions;
    Decoded_ASL = zeros(size(ASL_im,1),size(ASL_im,2),size(ASL_im,3),DecodedDataSize);
    
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
    
    % Save 4D nifti file
    xASL_io_SaveNifti(Encoded_ASL, Encoded_ASL, Decoded_ASL)
    
    % For visualization/testing -> remove this afterwards
%     figure 
%     image(Decoded_ASL(:,:,12,5));
    
end