% ------------------------------------------------------------------------------------------------------------------------------------------------------
% SHORT DESCRIPTION: Compiles all the MEX files
%
% FORMAT: xASL_compile_all
%
% INPUT:  n/a
%
% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This script can be used to compile all the MEX files.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     n/a
%
% REFERENCES:  n/a
% __________________________________
% Copyright © 2015-2018 ExploreASL
%
% 2018-12-18, Jan Petr
%

% Setup for C++ compiling

if ispc
    setenv('MW_MINGW64_LOC','c:\MinGW\64\')
end

mex -setup C++

mex xASL_mex_chamfers3D.c
mex xASL_mex_conv3Dsep.c