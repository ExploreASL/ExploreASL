function [Decoded_ASL] = xASL_im_HadamardDecoding(HadamardType, sec)
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


    % decides if Walsh or Natural Hadamard Matrix
    % for different sizes of matrices

    if sec == 1 %natural ordering
        
        if HadamardType == 4
            
            ASL_im = xASL_io_Nifti2Im('C:path to\ASL4D.nii.gz'); % Load time-series nifti
            dim4 = size(ASL_im, 4); %just to check the size dim4 (no. of volumes = 64)
            HadamardType = 4;
            nTE=8;
            rep=2;
            decoded_ASL_size = (HadamardType-1)*nTE*rep; %just to check size of decoded data
            
            %decoding scheme HAD-4
            %abs(a-b+c-d) = TI-1
            %abs(a-b-c+d) = TI-2
            %abs(a+b-c-d) = TI-3
            
            j=1;
            for x = 1:HadamardType:dim4
                
                Decoded_ASL(:,:,:,j) = ASL_im(:,:,:,(x)) - ASL_im(:,:,:,(x+1)) + ASL_im(:,:,:,(x+2)) - ASL_im(:,:,:,(x+3));
                Decoded_ASL(:,:,:,j+1) = ASL_im(:,:,:,(x)) - ASL_im(:,:,:,(x+1)) - ASL_im(:,:,:,(x+2)) + ASL_im(:,:,:,(x+3));
                Decoded_ASL(:,:,:,j+2) = ASL_im(:,:,:,(x)) + ASL_im(:,:,:,(x+1)) - ASL_im(:,:,:,(x+2)) - ASL_im(:,:,:,(x+3));
                
                j=j+(HadamardType-1);
            end
            
            figure %for visualization
            image(Decoded_ASL(:,:,12,1));
            
            % write the 4D nifti file
            Decoded_ASL_nii = make_nii(Decoded_ASL)
            save_nii(Decoded_ASL_nii, 'path to...\example_decoding_had_4.nii.gz')
            
        elseif HadamardType == 8
            
            %decoding scheme HAD-8
            %abs(a-b+c-d+e-f+g-h) = TI-1
            %abs(a-b+c-d-e+f-g+h) = TI-2
            %abs(a-b-c+d-e+f+g-h) = TI-3
            %abs(a-b-c+d+e-f-g+h) = TI-4
            %abs(a+b-c-d+e+f-g-h) = TI-5
            %abs(a+b-c-d-e-f+g+h) = TI-6
            %abs(a+b+c+d-e-f-g-h) = TI-7
            
            ASL_im = xASL_io_Nifti2Im('path to...\had_8_multite_encoded.nii.gz'); % Load time-series nifti
            dim4 = size(ASL_im, 4); %just to check the size dim4 (no. of volumes = 64)
            HadamardType = 8;
            nTE=8;
            rep=1;
            decoded_ASL_size = (HadamardType-1)*nTE*rep;
            
            j=1;
            for x = 1:HadamardType:dim4
                
                Decoded_ASL(:,:,:,j) = ASL_im(:,:,:,(x)) - ASL_im(:,:,:,(x+1)) + ASL_im(:,:,:,(x+2)) - ASL_im(:,:,:,(x+3)) + ASL_im(:,:,:,(x+4)) - ASL_im(:,:,:,(x+5)) + ASL_im(:,:,:,(x+6)) - ASL_im(:,:,:,(x+7));
                Decoded_ASL(:,:,:,j+1) = ASL_im(:,:,:,(x)) - ASL_im(:,:,:,(x+1)) + ASL_im(:,:,:,(x+2)) - ASL_im(:,:,:,(x+3)) - ASL_im(:,:,:,(x+4)) + ASL_im(:,:,:,(x+5)) - ASL_im(:,:,:,(x+6)) + ASL_im(:,:,:,(x+7));
                Decoded_ASL(:,:,:,j+2) = ASL_im(:,:,:,(x)) - ASL_im(:,:,:,(x+1)) - ASL_im(:,:,:,(x+2)) + ASL_im(:,:,:,(x+3)) - ASL_im(:,:,:,(x+4)) + ASL_im(:,:,:,(x+5)) + ASL_im(:,:,:,(x+6)) - ASL_im(:,:,:,(x+7));
                Decoded_ASL(:,:,:,j+3) = ASL_im(:,:,:,(x)) - ASL_im(:,:,:,(x+1)) - ASL_im(:,:,:,(x+2)) + ASL_im(:,:,:,(x+3)) + ASL_im(:,:,:,(x+4)) - ASL_im(:,:,:,(x+5)) - ASL_im(:,:,:,(x+6)) + ASL_im(:,:,:,(x+7));
                Decoded_ASL(:,:,:,j+4) = ASL_im(:,:,:,(x)) + ASL_im(:,:,:,(x+1)) - ASL_im(:,:,:,(x+2)) - ASL_im(:,:,:,(x+3)) + ASL_im(:,:,:,(x+4)) + ASL_im(:,:,:,(x+5)) - ASL_im(:,:,:,(x+6)) - ASL_im(:,:,:,(x+7));
                Decoded_ASL(:,:,:,j+5) = ASL_im(:,:,:,(x)) + ASL_im(:,:,:,(x+1)) - ASL_im(:,:,:,(x+2)) - ASL_im(:,:,:,(x+3)) - ASL_im(:,:,:,(x+4)) - ASL_im(:,:,:,(x+5)) + ASL_im(:,:,:,(x+6)) + ASL_im(:,:,:,(x+7));
                Decoded_ASL(:,:,:,j+6) = ASL_im(:,:,:,(x)) + ASL_im(:,:,:,(x+1)) + ASL_im(:,:,:,(x+2)) + ASL_im(:,:,:,(x+3)) - ASL_im(:,:,:,(x+4)) - ASL_im(:,:,:,(x+5)) - ASL_im(:,:,:,(x+6)) - ASL_im(:,:,:,(x+7));

                j=j+(HadamardType-1);
            end
            
            figure %for visualization
            image(Decoded_ASL(:,:,12,40));
            
            % write the 4D nifti file
            Decoded_ASL_nii = make_nii(Decoded_ASL);
            save_nii(Decoded_ASL_nii, 'path to...\example_decoding_had_8.nii.gz')
            
        else
            
            %output = decoding scheme for hadamard-12;
            
        end
        
    elseif asl == 2 % walsh ordered
        
        if HadamardType == 4
            
            %output = decoding scheme for hadamard-4;
            
        elseif HadamardType == 8
            
            %output = decoding scheme for hadamard-8;
            
        else
            
            %output = decoding scheme for hadamard-12;
            
        end
        
    end

end 

