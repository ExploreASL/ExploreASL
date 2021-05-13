function xASL_qc_PrintOrientation(niftiList, outputDir, outputFile)
% xASL_qc_PrintOrientation List NifTI orientation matrix
%
% FORMAT:       xASL_qc_PrintOrientation(niftiList, outputDir, outputFile);
% 
% INPUT:
% niftiList     - cell array with paths to NIfTI files. Can be .nii or
%                 .nii.gz. (REQUIRED)
% outputDir     - folder where the output TSV file will be saved (REQUIRED)
% outputFile    - filename as suffix to output TSV file (REQUIRED)
%
% OUTPUT:       n/a
% OUTPUT FILE:
%               - xASL_qc_PrintOrientation*.tsv containing several orientation
%               parameters
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function lists NifTI orientation matrices before and after
% image processing, respectively nii.mat0 and nii.mat. In ExploreASL this
% is used for QC to detect accidental left-right flips, as these can occur
% unnoticed as the brain structure appears relatively symmetrical.
% This can be detected by negative determinants.
% Also, this can be used to detect any significant differences in
% acquisition or image processing.
% 
% So orientation parameters and determinants should be similar across
% all scans from single scanner/coil, and registration should not
% give a relatively negative determinant.
% Results are saved in a TSV file
%
% This functions performs the following steps:
% 1. Print the header
% 2. Load the data
% 3. Print original orientation matrix
% 4. Print current orientation matrix
% 5. Print registration transformation matrix
% 6. Print FileName
% 7. Get statistics (mean & SD)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      xASL_qc_PrintOrientation(niftiList, outputDir, outputFile)
% __________________________________
% Copyright 2015-2021 ExploreASL

%% Administration TSV
savePath = fullfile(outputDir,['xASL_qc_PrintOrientation_' outputFile '.tsv']);

fclose all;
xASL_delete(savePath);

Summary_fid = fopen(savePath,'wt');

%% 1. Print the header
HeaderPrint = {'Transl_X' 'Transl_Y' 'Transl_Z' 'Rot_X' 'Rot_Y' 'Rot_Z' 'Scale_X' 'Scale_Y' 'Scale_Z' 'Shear_X' 'Shear_Y' 'Shear_Z'};
EndPrint = {'DetOrigin' 'DetCurrent' 'DetTransform'};

for iE=1:length(EndPrint)
    for iH=1:length(HeaderPrint)
        fprintf(Summary_fid,'%s\t',HeaderPrint{iH});
    end
    fprintf(Summary_fid,'%s\t',EndPrint{iE});
end
fprintf(Summary_fid,'%s','FileName');

%% 2. Load the data
if ~iscell(niftiList)
    niftiList = {niftiList}; % temporary fix for individual files rather than lists of files= {reg_exp_INPUT}; % temporary fix for individual files rather than lists of files
end

if ~isempty(niftiList)
    fprintf('%s', 'Summarizing orientation parameters:   ');
    
    for iNifti=1:numel(niftiList)
        xASL_TrackProgress(iNifti, numel(niftiList));
        tNII = xASL_io_ReadNifti(niftiList{iNifti});

        fprintf( Summary_fid, '\n');

        %% 3. Print original orientation matrix
        TempVec         = spm_imatrix(tNII.mat0);
        for iI=1:12;    fprintf(Summary_fid,'%s\t', num2str(round(100*TempVec(iI))/100) );end
        matrix2(:,iNifti)   = TempVec;
        det2(iNifti)        = det(tNII.mat0(1:3,1:3));
        fprintf(Summary_fid,'%s\t', num2str( det2(iNifti) ) );

        %% 4. Print current orientation matrix
        TempVec         = spm_imatrix(tNII.mat);
        for iI=1:12;    fprintf(Summary_fid,'%s\t', num2str(round(100*TempVec(iI))/100) );end
        matrix1(:,iNifti)   = TempVec;
        det1(iNifti)        = det(tNII.mat(1:3,1:3));
        fprintf(Summary_fid,'%s\t', num2str( det1(iNifti) ) );

        %% 5. Print registration transformation matrix
        TempReg         = tNII.mat/tNII.mat0;
        TempVec         = spm_imatrix(TempReg);
        for iI=1:12;    fprintf(Summary_fid,'%s\t', num2str(round(100*TempVec(iI))/100) );end
        matrix0(:,iNifti)   = TempVec;
        det0(iNifti)        = det(TempReg(1:3,1:3));
        fprintf(Summary_fid,'%s\t', num2str( det0(iNifti) ) );

        %% 6. Print FileName
        fprintf(Summary_fid,'%s\t', niftiList{iNifti} );

    end
    
    fprintf('\n%s\n', ['into: ' savePath]);

    %% ----------------------------------------------------------------------------
    %% 7. Get statistics (mean & SD)

    fprintf(Summary_fid,'\n\n');

    TotalOri = matrix2';
    TotalCur = matrix1';
    TotalTra = matrix0';


    % Print means
    for iI=1:12;    fprintf(Summary_fid,'%s\t', num2str(round(100*mean(TotalOri(:,iI)))/100) );end
    fprintf(Summary_fid,'%s\t', num2str(round(100*mean(det2)/100)) );

    for iI=1:12;    fprintf(Summary_fid,'%s\t', num2str(round(100*mean(TotalCur(:,iI)))/100) );end
    fprintf(Summary_fid,'%s\t', num2str(round(100*mean(det1)/100)) );

    for iI=1:12;    fprintf(Summary_fid,'%s\t', num2str(round(100*mean(TotalTra(:,iI)))/100) );end
    fprintf(Summary_fid,'%s\t', num2str(round(100*mean(det0)/100)) );

    fprintf(Summary_fid,'%s\t', 'Mean' );
    fprintf(Summary_fid,'\n');

    % Print SD
    for iI=1:12;    fprintf(Summary_fid,'%s\t', num2str(round(100*std(TotalOri(:,iI)))/100) );end
    fprintf(Summary_fid,'%s\t', num2str(round(100*std(det2)/100)) );

    for iI=1:12;    fprintf(Summary_fid,'%s\t', num2str(round(100*std(TotalCur(:,iI)))/100) );end
    fprintf(Summary_fid,'%s\t', num2str(round(100*std(det1)/100)) );

    for iI=1:12;    fprintf(Summary_fid,'%s\t', num2str(round(100*std(TotalTra(:,iI)))/100) );end
    fprintf(Summary_fid,'%s\t', num2str(round(100*std(det0)/100)) );

    fprintf(Summary_fid,'%s\t', 'StDev' );
    fprintf(Summary_fid,'\n');


end
fclose all;


end