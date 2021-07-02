function [sec] = xASL_im_HadamardCaseByCase(hadamard_type)
%xASL_im_HadamardCaseByCase Case by case function for determining hadamard type and respective decoding scheme
%
% FORMAT:       [sec] = xASL_im_HadamardCaseByCase(hadamard_type)
% 
% INPUT:        hadamard_type - ...
%
% OUTPUT:       sec           - ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Case by case function for determining hadamard type and respective decoding scheme.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      n/a
% __________________________________
% Copyright 2015-2021 ExploreASL

    % Important feedback from Beatriz:
    % Add warning if it is a string or not and also if it is emplty or not

    if ( hadamard_type == natural)
        sec = 1;
    else
        sec = 2;
    end

end
