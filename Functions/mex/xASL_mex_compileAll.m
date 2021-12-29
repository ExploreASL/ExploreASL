function xASL_mex_compileAll()
%xASL_mex_compileAll Compiles the whole xASL mex code
%
% FORMAT: xASL_compile_all
%
% INPUT:  n/a
%
% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This script can be used to compile all the MEX files that are in the mex directory, in the SPMmodified
%              and also the DCMTK files.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     n/a
%
% REFERENCES:  n/a
% __________________________________
% Copyright © 2015-2021 ExploreASL

% Obtain the current path
pathCurrent = matlab.desktop.editor.getActiveFilename;

% Remove the filename
[pathCurrent,~,~] = fileparts(pathCurrent);

% Go to the basic ExploreASL folder
pathXASL = pathCurrent(1:end-4);

% Setup for C++ compiling	
if ispc
    setenv('MW_MINGW64_LOC','c:\MinGW\64\')
end

mex -setup C++

% Compile the xASL files in the SPM folder
cd(fullfile(pathXASL,'External','SPMmodified','xASL'));
mex xASL_mex_chamfers3D.c
mex xASL_mex_conv3Dsep.c

% Compile the xASL files in the mex folder
cd(fullfile(pathXASL,'mex'));
mex xASL_mex_conv3DsepGauss.c
mex xASL_mex_dilate_erode_3D.c
mex xASL_mex_dilate_erode_single.c
mex xASL_mex_JointHist.c

% Compile the xASL DICOM reading using the DCMTK library
cd(fullfile(pathXASL,'External','DCMTK'));
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

% Compute the special SPM functions we use and modify
cd(fullfile(pathXASL,'External','SPMmodified','src'))
mex -O -largeArrayDims spm_jsonread.c jsmn.c -DJSMN_PARENT_LINKS

end
