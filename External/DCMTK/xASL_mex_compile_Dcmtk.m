% ------------------------------------------------------------------------------------------------------------------------------------------------------
% SHORT DESCRIPTION: Compiles the DCMTK MEX files
%
% FORMAT: xASL_compile_all
%
% INPUT:
% OUTPUT:
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: 
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:
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

% Compiles the DICOM reading MEX - calls with a different make-file for each system
if ismac
	mex -v -f xASL_mex_DcmtkMac.xml xASL_mex_DcmtkRead.cpp
else
    if isunix
        mex -v -f xASL_mex_DcmtkUnix.xml xASL_mex_DcmtkRead.cpp
    end
    
    if ispc
	mex -v -f xASL_mex_DcmtkWin.xml xASL_mex_DcmtkRead.cpp
    end
end