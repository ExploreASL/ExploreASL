function xASL_mex_compileAll()
%xASL_mex_compileAll Compiles all xASL-related c-code into MEX
%
% FORMAT: xASL_compile_all
%
% INPUT:  n/a
%
% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This script can be used to compile the following MEX files in the:
%              1. SPM xASL subfolder
%              2. xASL proper
%              3. DCMtK
%              4. SPM proper
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     xASL_compile_all;
% __________________________________
% Copyright © 2015-2023 ExploreASL


%% Admin
% Obtain the current path
pathCurrent = matlab.desktop.editor.getActiveFilename;

% Remove the filename
[pathCurrent,~,~] = fileparts(pathCurrent);

% Go to the basic ExploreASL folder
pathXASL = pathCurrent(1:end-14);

% Detect environment
if ispc
    fprintf('%s\n', 'Windows PC environment detected, trying recompiling...');
    if ~exist('c:\MinGW\64\', 'dir')
        error('No MinGW compilation tool found at the usual location, skipping');
    end
   
    setenv('MW_MINGW64_LOC', 'c:\MinGW\64\');
    environmentIs = 'w64';
elseif ismac
    [~, result] = system('uname -m');
    if strcmp(strtrim(result),'arm64')
        fprintf('%s\n', 'mac arm64 (Apple silicon M1/M2/etc processor) detected, trying recompiling...');
        environmentIs = 'maca64';
    else
        fprintf('%s\n', 'mac i64 (intel processor) detected, trying recompiling...');
        environmentIs = 'maci64';
    end
elseif isunix
    fprintf('%s\n', 'Linux detected, trying recompiling...');
    environmentIs = 'a64';
else
    error('Unknown or not supported environment, stopping compiling');
end

% Setup for C++ compiling
mex -setup C++


%% 1. Compile the xASL files in the SPM folder
fprintf('%s\n', 'Compiling the SPM xASL folder...');
cd(fullfile(pathXASL,'External','SPMmodified','xASL'));
mex xASL_mex_chamfers3D.c
mex xASL_mex_conv3Dsep.c


%% 2. Compile the xASL files in the mex folder
fprintf('%s\n', 'Compiling the xASL proper folder...');
cd(fullfile(pathXASL,'Functions','mex'));
mex xASL_mex_conv3DsepGauss.c
mex xASL_mex_dilate_erode_3D.c
mex xASL_mex_dilate_erode_single.c
mex xASL_mex_JointHist.c


%% 3. Compile the xASL DICOM reading using the DCMTK library
fprintf('%s\n', 'Compiling the DCMtK folder...');
cd(fullfile(pathXASL,'External','DCMTK'));

switch environmentIs
    case 'maca64'
        mex -v -f xASL_mex_DcmtkMacARM64.xml xASL_mex_DcmtkRead.cpp
    case 'maci64'
        mex -v -f xASL_mex_DcmtkMacINTEL64.xml xASL_mex_DcmtkRead.cpp
    case 'a64'
        mex -v -f xASL_mex_DcmtkUnix.xml xASL_mex_DcmtkRead.cpp
    case 'w64'
	    mex -v -f xASL_mex_DcmtkWin.xml xASL_mex_DcmtkRead.cpp
end

end
