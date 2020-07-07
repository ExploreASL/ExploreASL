function xASL_qc_PrintOrientation(DIR, reg_exp_INPUT,OUTPUT_DIR,Name)
% xASL_qc_PrintOrientation Check orientation of niftis, useful to detect
% accidental left-right flips (all other flips will be visible).
% translations, rotations or shears are not to be worried about,
% only negative zooms. This can be detected by negative determinants.
% So orientation parameters and determinants should be similar across
% all scans from single scanner/coil, and registration should not
% give negative determinant.
%
% CAVE: will search recursively through directories for niftis that fullfill reg_exp_INPUT!
%
% FORMAT:       xASL_qc_PrintOrientation(DIR, reg_exp_INPUT,OUTPUT_DIR,Name);
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Check orientation of niftis, useful to detect
%               accidental left-right flips (all other flips will be visible).
%               translations, rotations or shears are not to be worried about,
%               only negative zooms. This can be detected by negative determinants.
%               So orientation parameters and determinants should be similar across
%               all scans from single scanner/coil, and registration should not
%               give negative determinant.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL

%% Administration CSV
fclose all;
xASL_adm_DeleteFileList(OUTPUT_DIR, ['xASL_qc_PrintOrientation_' Name '\.(csv|tsv)'],[],[0 Inf]);
SaveFile = fullfile(OUTPUT_DIR,['xASL_qc_PrintOrientation_' Name '.tsv']);
Summary_fid = fopen(SaveFile,'wt');

%% Header
HeaderPrint = {'Transl_X' 'Transl_Y' 'Transl_Z' 'Rot_X' 'Rot_Y' 'Rot_Z' 'Scale_X' 'Scale_Y' 'Scale_Z' 'Shear_X' 'Shear_Y' 'Shear_Z'};
EndPrint = {'DetOrigin' 'DetCurrent' 'DetTransform'};

for iE=1:length(EndPrint)
    for iH=1:length(HeaderPrint)
        fprintf(Summary_fid,'%s\t',HeaderPrint{iH});
    end
    fprintf(Summary_fid,'%s\t',EndPrint{iE});
end
fprintf(Summary_fid,'%s','FileName');

%% Get data
% FList = xASL_adm_GetFileList(DIR, reg_exp_INPUT, 'FPListRec', [0 Inf]);

FList = {reg_exp_INPUT}; % temporary fix for individual files rather than lists of files

if ~isempty(FList)

    for iF=1:length(FList)
        tNII = xASL_io_ReadNifti(FList{iF});

        fprintf( Summary_fid, '\n');

        % Print original orientation matrix
        TempVec         = spm_imatrix(tNII.mat0);
        for iI=1:12;    fprintf(Summary_fid,'%s\t', num2str(round(100*TempVec(iI))/100) );end
        matrix2(:,iF)   = TempVec;
        det2(iF)        = det(tNII.mat0(1:3,1:3));
        fprintf(Summary_fid,'%s\t', num2str( det2(iF) ) );

        % Print current orientation matrix
        TempVec         = spm_imatrix(tNII.mat);
        for iI=1:12;    fprintf(Summary_fid,'%s\t', num2str(round(100*TempVec(iI))/100) );end
        matrix1(:,iF)   = TempVec;
        det1(iF)        = det(tNII.mat(1:3,1:3));
        fprintf(Summary_fid,'%s\t', num2str( det1(iF) ) );

        % Print registration transformation matrix
        TempReg         = tNII.mat/tNII.mat0;
        TempVec         = spm_imatrix(TempReg);
        for iI=1:12;    fprintf(Summary_fid,'%s\t', num2str(round(100*TempVec(iI))/100) );end
        matrix0(:,iF)   = TempVec;
        det0(iF)        = det(TempReg(1:3,1:3));
        fprintf(Summary_fid,'%s\t', num2str( det0(iF) ) );

        % Print FileName
        fprintf(Summary_fid,'%s\t', FList{iF} );

    end

    %% ----------------------------------------------------------------------------
    %% Get statistics

    fprintf(Summary_fid,'\n\n');

    TotalOri        = matrix2';
    TotalCur        = matrix1';
    TotalTra        = matrix0';


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
