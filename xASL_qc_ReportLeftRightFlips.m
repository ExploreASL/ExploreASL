% Check for left-right flips for all NIfTIs inside a folder
dirRoot = '/Users/henk/surfdrive/HolidayPics/CICERO_Nolan/analysis';

niftiList = xASL_adm_GetFileList(dirRoot, '^.*\.nii$', 'FPListRec');
niftiList = niftiList(1:30);
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
        
