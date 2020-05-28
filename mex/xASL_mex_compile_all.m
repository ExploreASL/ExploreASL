% ------------------------------------------------------------------------------------------------------------------------------------------------------
% SHORT DESCRIPTION: Compiles all the MEX files
%
% FORMAT: xASL_compile_all
%
% INPUT:
% OUTPUT:
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function compiles all MEX files inside this folder.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:
% __________________________________
% Copyright © 2015-2020 ExploreASL
%

% Setup for C++ compiling

if ispc
    setenv('MW_MINGW64_LOC','c:\MinGW\64\')
end

mex -setup C++

mex xASL_mex_conv3DsepGauss.c
mex xASL_mex_dilate_erode_3D.c
mex xASL_mex_dilate_erode_single.c
mex xASL_mex_JointHist.c