function [identical,RMSE,minDiff,maxDiff,dimCheck] = xASL_im_CompareNiftis(pathA,pathB,bVerbose)
%xASL_im_CompareNiftis Compare two niftis. Untouched comparison based on copies.
%
% FORMAT: [identical,RMSE] = xASL_im_CompareNiftis(pathA,pathB)
%
% INPUT:
%        pathA         - path nifti A (CHAR ARRAY, REQUIRED)
%        pathB         - path nifti B (CHAR ARRAY, REQUIRED)
%        bVerbose      - print information (BOOLEAN, OPTIONAL, DEFAULT=true)
%
% OUTPUT:
%        identical     - True if both nifti images are exactly the same
%        RMSE          - Root mean square error
%        minDiff       - abs(min(imageA) - min(imageB))
%        maxDiff       - abs(max(imageA) - max(imageB))
%        dimCheck      - True if both images have the same dimensions
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% DESCRIPTION:      Compare two niftis. Untouched comparison based on copies.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          [identical,RMSE] = xASL_im_CompareNiftis(pathA,pathB);
%
% __________________________________
% Copyright @ 2015-2021 ExploreASL


    %% Input check
    
    % Defaults
    identical = false;
    RMSE = inf;
    minDiff = inf;
    maxDiff = inf;
    dimCheck = false;
    
    % Arguments
    if nargin<2
        error('Image paths missing...');
    end
    if nargin<3
        bVerbose = true;
    end
    
    % Check existence
    if ~xASL_exist(pathA) || ~xASL_exist(pathB)
        error('Image A or image B does not exist...');
    end
    
    % Get file names
    [baseA,nameNiftiA,extensionA] = xASL_fileparts(char(pathA));
    [baseB,nameNiftiB,extensionB] = xASL_fileparts(char(pathB));
    
    % Check extensions
    if (~strcmp(extensionA,'.nii.gz') && ~strcmp(extensionA,'.nii')) ...
    || (~strcmp(extensionB,'.nii.gz') && ~strcmp(extensionB,'.nii'))
        error('Image A or image B have a wrong file type...');
    end
    
    % Make sure pathA is not equal to pathB
    if strcmp(pathA,pathB)
        error('You inserted the same image path twice...');
    end
    
    % Make sure the file size is not zero (corrupted files etc.)
    structA = dir(pathA);
    structB = dir(pathB);
    sizeA = structA.bytes;
    sizeB = structB.bytes;
    if ~(sizeA>1) || ~(sizeB>1)
        warning('The niftis seem damanged. The file size is too small...');
        return
    end
    
    %% Untouched loading of niftis (copy, then load so that original ones arent unzipped, then delete copies)
    
    % Define copy names
    nameCopyA = [nameNiftiA '_tmp_copy' extensionA];
    nameCopyB = [nameNiftiB '_tmp_copy' extensionB];
    pathCopyA = fullfile(baseA, nameCopyA);
    pathCopyB = fullfile(baseB, nameCopyB);
    
    % Copy files
    xASL_Copy(char(pathA),pathCopyA);
    xASL_Copy(char(pathB),pathCopyB);
    
    % Load image copies
    imageA = xASL_io_Nifti2Im(pathCopyA);
    imageB = xASL_io_Nifti2Im(pathCopyB);
    
    % Delete copies
    xASL_delete(pathCopyA);
    xASL_delete(pathCopyB);
    
    %% Compare the niftis
    
    % Identical check
    if isequal(imageA,imageB)
        identical = true;
    end
    
    % Get RMSE
    if isequal(size(imageA),size(imageB))
        dimCheck = true;
        RMSE = sqrt(mean((imageA(:) - imageB(:)).^2))*2/sqrt(mean(abs(imageA(:)) + abs(imageB(:))).^2);
    end
    
    % Min/Max A
    minA = min(imageA(:));
    maxA = max(imageA(:));
    
    % Min/Max B
    minB = min(imageB(:));
    maxB = max(imageB(:));
    
    minDiff = abs(minA-minB);
    maxDiff = abs(maxA-maxB);
    
    %% Print statements
    if bVerbose
        if identical
            fprintf('Identical:       true\n');
        else
            fprintf('Identical:       false\n');
        end
        fprintf('RMSE:            %.2f\n', RMSE);
        fprintf('Min. difference: %.2f\n', minDiff);
        fprintf('Max. difference: %.2f\n', maxDiff);
    end

end




