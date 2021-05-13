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
% DESCRIPTION: This submodule registers T1w linearly to the center of MNI space, a.k.a. ACPC alignment
% The same transformation is applied to all other related scans (ASL4D, M0, FLAIR, etc.)
% This facilitates MNI-based algorithms (e.g. SPM-based segmentation), and allows for visual QC with all images
% roughly in the same space. This submodule first clips high values that can bias the registration algorithm, then performs
% a center of mass-based ACPC alignment, and then several iterations of SPM coregistration.
% Assuming that this submodule is run at the start of ExploreASL, all NIfTI orientation matrices are restored before running the registration.
%
% All NifTIs are found recursively (i.e. in the folder and its subfolders),
% irregardless of them being .nii or .nii.gz.

% EXAMPLE: xASL_qc_ReportLeftRightFlips('/Users/henk/surfdrive/HolidayPics/CICERO_Nolan/analysis');
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2019 ExploreASL


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

    if ~isempty(LR_flip_YesNo)
        fprintf('%s\n', 'Left-right flip(s) detected in:');

        for iFlip=1:numel(LR_flip_YesNo)
            fprintf('%s\n', niftiList{LR_flip_YesNo(iFlip)});
        end
        fprintf('\n');
    end
    
    if bZip
        xASL_adm_GzipAllFiles(dirRoot);
    end
        
end