function [sec] = xASL_im_HadamardCaseByCase(hadamard_matrix)
%xASL_im_HadamardCaseByCase Case by case function for determining hadamard type and respective decoding scheme
%
% FORMAT:       [sec] = xASL_im_HadamardCaseByCase(hadamard_type)
% 
% INPUT:        hadamard_matrix - should be a string, either 'Natural' or % 'Walsh'
%
% OUTPUT:       sec           - defines the section to continue, either % 1(for Natural) or 2(for Walsh)
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Case by case function for determining hadamard type and respective decoding scheme.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      n/a
% __________________________________
% Copyright 2015-2021 ExploreASL

  
    if ~ischar(hadamard_matrix) || isempty(hadamard_matrix) 
        warning('Hadamard matrix is not a string or is empty')
        
    elseif ~strcmpi(hadamard_matrix,'natural') || ~strcmpi(hadamard_matrix,'walsh')
        warning('Hadamard matrix is not defined')
    
    elseif %HadamatrixMatrix is correctlydefined
        if strcmpi(hadamard_matrix,'natural')
            sec = 1;
        else
            sec = 2;
        end
    end
     

   

end
