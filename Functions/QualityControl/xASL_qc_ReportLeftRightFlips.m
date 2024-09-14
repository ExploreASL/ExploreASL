function xASL_qc_ReportLeftRightFlips(dirRoot, bZip)
%xASL_qc_ReportLeftRightFlips Identify & report LR-flips for all NIfTIs inside a folder
%
% FORMAT: xASL_qc_ReportLeftRightFlips(dirRoot [, bZip])
%
% INPUT:
%   dirRoot     - folder containing NIfTIs (REQUIRED)
%   bzip        - true for zipping NIfTIs after they were unzipped
%                 (OPTIONAL, DEFAULT = true)
%
% OUTPUT:
% Prints a list on the screen with NIfTIs that contain left-right flips
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function identifies and reports illegal left-right
% flips for image matrices within a NIfTI. This can be useful as these are
% not readily observed by the human eye, as the left and right hemispheres
% are too symmetrical by default.
%
% All NifTIs are found recursively (i.e. in the folder and its subfolders),
% irregardless of them being .nii or .nii.gz.
%
% EXAMPLE: xASL_qc_ReportLeftRightFlips('/Users/henk/surfdrive/HolidayPics/CICERO_Nolan/analysis');
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2021 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

    if nargin<1 || isempty(dirRoot)
        error('Root folder missing');
    elseif nargin<2 || isempty(bZip)
        bZip = true;
    end
    niftiList = xASL_adm_GetFileList(dirRoot, '^.*\.nii$', 'FPListRec');
    fileName = 'AllNiftis';
    xASL_qc_PrintOrientation(niftiList, dirRoot, fileName);
    PathOrientationResults = fullfile(dirRoot, ['xASL_qc_PrintOrientation_' fileName '.tsv']);
    LR_flip_YesNo = xASL_im_DetermineFlip(PathOrientationResults);
    savePath = fullfile(dirRoot, 'ReportLeftRightFlips.tsv');
    
    if ~isempty(LR_flip_YesNo)
        
        fclose('all');
        xASL_delete(savePath);
        fID = fopen(savePath, 'wt');
        
        fprintf('%s\n', ['Left-right flip(s) detected: ' savePath]);
        
        for iFlip=1:numel(LR_flip_YesNo)
            fprintf(fID, '%s\n', niftiList{LR_flip_YesNo(iFlip)});
        end
        fprintf('\n');
        
        fclose('all');
    end
    
    if bZip
        xASL_adm_GzipAllFiles(dirRoot);
    end
        
end