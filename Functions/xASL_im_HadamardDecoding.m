function [Decoded_ASL] = xASL_im_HadamardDecoding(TimeEncodedMatrixType, TimeEncodedMatrixSize, DecodingMatrix_input)
%xASL_im_HadamardDecoding Hadamard-4 & Hadamard-8 Decoding
%
% FORMAT:       [Decoded_ASL] = xASL_im_HadamardDecoding(HadamardType, sec)
%
% INPUT:        HadamardType  - 4 for Hadamard4
%                             - 8 for Hadamard8
%                             - 12 for Hadamard12
%               sec           - output from xASL_im_HadamardCaseByCase
%
% OUTPUT:       Decoded_ASL   - ...
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Hadamard-4 & Hadamard-8 Decoding.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      n/a
% __________________________________
% Copyright 2015-2021 ExploreASL


% #### For Hadamard Decoding Matrix ####

if strcmp(TimeEncodedMatrixType,'Hadamard')
    
    if TimeEncodedMatrixSize == 4
        
        ASL_im = xASL_io_Nifti2Im('C:\Users\amahroo\Documents\GitHub\Siemens_FME_Hadamard\Siemens_FME_Hadamard\derivatives\ExploreASL\sub-Sub1\ASL_1\ASL4D.nii'); % Load time-series nifti
        
        if ~isempty(DecodingMatrix_input) %if there's a Decoding Matrix from the data
            DecodingMatrix = DecodingMatrix_input;
        else
            DecodingMatrix = [1 -1 1 -1;
                1 -1 -1 1;
                1 1  -1 -1];
        end
        
    elseif TimeEncodedMatrixSize == 8
        
        ASL_im = xASL_io_Nifti2Im('C:\Users\amahroo\Desktop\multite_study\ExploreASL\multite_had_8\had_8_multite_encoded.nii.gz'); % Load time-series nifti
        
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

    nTE = 8;                                         %to be replaced by x.HadamardNumberTE
    nDecodedTI = HadamardSize-1;                     %to be replaced by x.HadamardMatrixSize - 1
    nDecodedVol = nTE * nDecodedTI;
    nEncodedVol = size(ASL_im, 4);                   %calculating no. of repeats - total volumes of encoded data
    EncodedDataSize = HadamardSize * nTE;            %expected data size
    nRep = nEncodedVol / EncodedDataSize;
    
    DecodedDataSize = nDecodedVol * nRep;
    Decoded_ASL = zeros(size(ASL_im,1),size(ASL_im,2),size(ASL_im,3),DecodedDataSize);
    
    idx=0;
    for repetition = 1:nRep
        for TE = 1:nTE
            for TI = 1:nDecodedTI
                
                indexPositive = find(DecodingMatrix(TI,:)==1);
                indexNegative = find(DecodingMatrix(TI,:)==-1);
                Decoded_ASL(:,:,:,((TE-1)*nDecodedTI+TI)+(repetition-1)*nDecodedVol) = mean(ASL_im(:,:,:,(indexPositive+idx)),4) - mean(ASL_im(:,:,:,(indexNegative+idx)),4);
                
            end
            idx=idx+HadamardSize;
        end
    end
    
    %write 4D nifti file
    Decoded_ASL_nii = make_nii(Decoded_ASL);
    save_nii(Decoded_ASL_nii, 'C:\Users\amahroo\Desktop\multite_study\ExploreASL\decoding_test\example_decoding.nii.gz');
    
    figure %for visualization
    image(Decoded_ASL(:,:,12,5));
    
end