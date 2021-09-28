function [pathOut] = xASL_spm_admin(pathIn, bPadComma1)
%xASL_spm_admin Manage SPM input path format
%
% FORMAT: xASL_spm_admin(pathIn[, bPadComma1])
%
% INPUT:
%   pathIn     - path to NifTI image, including full foldername and extension (REQUIRED)
%   bPadComma1 - boolean specifying if the ',1' is added, 
%                to assign the first image of an SPM 4D array (OPTIONAL, DEFAULT = FALSE)
%
% OUTPUT:
%   pathOut    - corrected path to NifTI image
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This SPM admin function takes a NIfTI path and does a few
% checks to make this valid to SPM. It accepts both .nii and .nii.gz.
% It runs the following steps:
%
% 1. Unzip .nii.gz
% 2. Convert char to cell
% 3. Add ',1' suffix
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_spm_admin('/MyStudy/Subject1/T1.nii.gz');
% __________________________________
% Copyright 2015-2021 ExploreASL


    if nargin<1 || isempty(pathIn)
        error('Too few input arguments');
    elseif ~ischar(pathIn)
        error('pathIn should be a char array');
    end    

    if nargin<2 || isempty(bPadComma1)
        bPadComma1 = true; % default
    end

    %% -----------------------------------------------------
    %% 1. Unzip .nii.gz
    [~, pathIn] = xASL_io_ReadNifti(pathIn);

    
    %% -----------------------------------------------------
    %% 2. Convert char to cell
	if ~iscell(pathIn)
		tempL = pathIn;
		IMnameOut{1} = tempL;
	else
		IMnameOut = pathIn;
    end

    pathIn = IMnameOut;

    
    %% -----------------------------------------------------
    %% 3. Add ',1' suffix
    if bPadComma1
        String2Pad = ',1';
    else
        String2Pad = '';
    end
    
    FoundStr = strfind(pathIn{1}, ',1');
    if isempty(FoundStr)
        pathIn{1} = [pathIn{1} String2Pad];
    else
        FoundStr = FoundStr(length(FoundStr));
        if FoundStr+1~=length(pathIn{1})
            pathIn{1} = [pathIn{1} String2Pad];
        end
    end

    pathOut = pathIn;
    
end
