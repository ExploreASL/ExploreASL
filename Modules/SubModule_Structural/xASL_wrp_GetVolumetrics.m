function xASL_wrp_GetVolumetrics(x)
%xASL_wrp_GetVolumetrics Submodule of ExploreASL Structural Module, that obtains volumes from the tissue segmentations
% (& FLAIR WMH segmentations if they exist)
%
% FORMAT: xASL_wrp_GetVolumetrics(x)
%
% INPUT:
%   x 	    - structure containing fields with all information required to run this submodule (REQUIRED)
%   x.P     - paths with NIfTIs for which this function should be applied to (REQUIRED)
%   x.WMHsegmAlg - Choose the LST algorithm (OPTIONAL 'LGA', DEFAULT = 'LPA')
%
% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This submodule computes the total volume for each of the tissue classes & stores them in a TSV file (per BIDS).
% This is computed from the native space segmentation derivatives (GM, WM & CSF), from which the ICV & relative volumes can be
% calculated. This is performed for CAT12 or SPM12 (whichever was used), and optionally for a WMH_SEGM.
%
% EXAMPLE: xASL_wrp_GetVolumetrics(x);
% __________________________________
% Copyright 2015-2020 ExploreASL


%% ------------------------------------------------------------------------------------------------
%% 1) File management & store tissue volumes
catVolFile = fullfile(x.D.TissueVolumeDir, ['cat_' x.P.STRUCT '_' x.P.SubjectID '.mat']); % CAT12 results
MatFile = fullfile(x.dir.SUBJECTDIR, [x.P.STRUCT '_seg8.mat']); % SPM12 results
SaveFile = fullfile(x.D.TissueVolumeDir, ['TissueVolume_' x.P.SubjectID '.csv']); % ExploreASL results
xASL_adm_CreateDir(x.D.TissueVolumeDir);

if xASL_exist(catVolFile, 'file') % for CAT12 segmentation
    catVol = load(catVolFile);
    if isfield(catVol.S,'subjectmeasures')
        GMvol = catVol.S.subjectmeasures.vol_abs_CGW(2)/1000; % volume was in cm^3, /1000 = dm^3 or Liter
        WMvol = catVol.S.subjectmeasures.vol_abs_CGW(3)/1000;
        CSFvol = catVol.S.subjectmeasures.vol_abs_CGW(1)/1000;
    else
        warning('catVol.S does not contain the subjectmeasures field...');
        GMvol = nan; WMvol = nan; CSFvol = nan;
    end

    FileID = fopen(SaveFile, 'wt');
    fprintf(FileID,'%s\n', 'File,GM_volume_L,WM_volume_L,CSF_volume_L');
    fprintf(FileID,'%s', [catVolFile ',' num2str(GMvol) ',' num2str(WMvol) ',' num2str(CSFvol)]);
    fclose(FileID);

elseif exist(MatFile, 'file') % for SPM12 segmentation

    % Make sure that filename is correct, otherwise this crashes.
    % This happens if the files have been run previously in different folder structure
    clear tpm image Affine lkp MT Twarp Tbias mg mn vr wp ll volumes
    load(MatFile, '-mat');
    for ii=1:6
        tpm(ii).fname = fullfile(x.D.SPMDIR, 'tpm', 'TPM.nii');
    end
    image.fname = fullfile(x.dir.SUBJECTDIR, [x.P.STRUCT '.nii']);

    save(MatFile, 'image', 'tpm', 'Affine', 'lkp', 'MT', 'Twarp', 'Tbias', 'mg', 'mn', 'vr', 'wp', 'll');

    matlabbatch{1}.spm.util.tvol.matfiles = {MatFile};
    matlabbatch{1}.spm.util.tvol.tmax = 6;
    matlabbatch{1}.spm.util.tvol.mask = {fullfile(x.D.SPMDIR, 'tpm', 'mask_ICV.nii,1')};
    matlabbatch{1}.spm.util.tvol.outf = SaveFile;

    xASL_adm_UnzipNifti(image.fname, 1); % unzip the NIfTI if needed
    spm_jobman('run', matlabbatch); % run the volumetrics
    
    % Edit the header names
    CSV = xASL_csvRead(SaveFile);
    CSV(1, 2:end) = {'GM_volume_L' 'WM_volume_L' 'CSF_volume_L' 'Tissue_volume_L' 'Bone_volume_L' 'Air_volume_L'};
    xASL_tsvWrite(CSV, [SaveFile(1:end-4) '.tsv'], 1);
    xASL_delete(SaveFile);
else
    warning('No CAT12/SPM12 volumetric result file found, need to repeat segmentation');
end

if exist(SaveFile, 'file')
    xASL_bids_csv2tsvReadWrite(SaveFile);
end


%% 2) File management & store volumes for WMH segmentation (if exists)
if xASL_exist(x.P.Path_WMH_SEGM, 'file')

    xASL_io_ReadNifti(x.P.Path_WMH_SEGM); % make sure it is unzipped

    % This function prints the WMH volume (mL) and number of WMH in a TSV file
    bVerbose = 1; % verbal output
    BinThresh = 0.05; % default setting = 0.5. we did a cleanup procedure already, so can allow a very low threshold
    MinimalLesionVolume = 0.015; % default (mL), determines the minimal size of each detected lesion
    
    % Delete any previous volume files:
    xASL_adm_DeleteFileList(x.dir.SUBJECTDIR, '^LST_tlv_.*\.(csv|tsv)$', [], [0 Inf]);
    xASL_adm_DeleteFileList(x.D.TissueVolumeDir, '^LST_tlv_.*\.(csv|tsv)$', [], [0 Inf]);
    xASL_adm_DeleteFileList(x.D.TissueVolumeDir, ['WMH_LST_(LGA|LPA)_' x.P.SubjectID '.(csv|tsv)'], [], [0 Inf]);
    
    % run the LST lesion volume calculator
    CurrDir = pwd;
    cd(x.D.TissueVolumeDir);
    [~, FileName] = ps_LST_tlv(x.P.Path_WMH_SEGM, ~bVerbose, BinThresh, MinimalLesionVolume);
    cd(CurrDir);
    % get current filename:
    FileName = fullfile(x.D.TissueVolumeDir, FileName);
    % define new filename:
    OutputPath = fullfile(x.D.TissueVolumeDir, ['WMH_LST_' x.WMHsegmAlg '_' x.P.SubjectID '.csv']);
    % Then move our file
    if ~exist(FileName, 'file')
         warning(['Missing:' FileName]);
    else
        xASL_Move(FileName, OutputPath, true);
        xASL_bids_csv2tsvReadWrite(OutputPath, true); % convert to tsv per BIDS
    end
end


end