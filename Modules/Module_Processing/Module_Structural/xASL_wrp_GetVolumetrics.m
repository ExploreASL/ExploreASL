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
% Copyright 2015-2023 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

%% ------------------------------------------------------------------------------------------------
%% 1) File management & store tissue volumes
fileCAT12Volume = ['cat_' x.P.STRUCT '_' x.P.SubjectID '.mat']; % CAT12 results
pathCAT12Volume = fullfile(x.D.TissueVolumeDir, fileCAT12Volume);
fileSPM12Volume = [x.P.STRUCT '_seg8.mat']; % SPM12 results
pathSPM12Volume = fullfile(x.dir.SUBJECTDIR, fileSPM12Volume);
SaveFileTSV = fullfile(x.D.TissueVolumeDir, ['TissueVolume_' x.P.SubjectID '.tsv']); % ExploreASL results
SaveFileCSV = fullfile(x.D.TissueVolumeDir, ['TissueVolume_' x.P.SubjectID '.csv']); % ExploreASL results
xASL_adm_CreateDir(x.D.TissueVolumeDir); % only create folder if it doesn't exist yet
%% For CAT12 segmentation
if xASL_exist(pathCAT12Volume, 'file')
    volumeCAT12 = load(pathCAT12Volume, '-mat');
    tableCAT12 = {'File' 'GM_volume_L' 'WM_volume_L' 'CSF_volume_L'};
    tableCAT12{2,1} = x.P.SubjectID;
    if isfield(volumeCAT12.S,'subjectmeasures')
        tableCAT12{2,2} = xASL_num2str(volumeCAT12.S.subjectmeasures.vol_abs_CGW(2)/1000); % GM volume volume was in cm^3, /1000 = dm^3 or Liter
        tableCAT12{2,3} = xASL_num2str(volumeCAT12.S.subjectmeasures.vol_abs_CGW(3)/1000); % WM volume
        tableCAT12{2,4} = xASL_num2str(volumeCAT12.S.subjectmeasures.vol_abs_CGW(1)/1000); % CSF volume
    else
        warning('catVol.S does not contain the subjectmeasures field...');
        tableCAT12{2,2:4} = 'n/a';
    end
    xASL_tsvWrite(tableCAT12, SaveFileTSV, 1);
%% For SPM12 segmentation
elseif xASL_exist(pathSPM12Volume, 'file')
    % Make sure that filename is correct, otherwise this crashes.
    % This happens if the files have been run previously in different folder structure
    % clear tpm image Affine lkp MT Twarp Tbias mg mn vr wp ll volumes
    
    load(pathSPM12Volume, '-mat');
    for ii=1:6
        tpm(ii).fname = fullfile(x.D.SPMDIR, 'tpm', 'TPM.nii');
    end
    image.fname = fullfile(x.dir.SUBJECTDIR, [x.P.STRUCT '.nii']);
    save(pathSPM12Volume, 'image', 'tpm', 'Affine', 'lkp', 'MT', 'Twarp', 'Tbias', 'mg', 'mn', 'vr', 'wp', 'll');
    matlabbatch{1}.spm.util.tvol.matfiles = {pathSPM12Volume};
    matlabbatch{1}.spm.util.tvol.tmax = 6;
    matlabbatch{1}.spm.util.tvol.mask = {fullfile(x.D.SPMDIR, 'tpm', 'mask_ICV.nii,1')};
    matlabbatch{1}.spm.util.tvol.outf = SaveFileCSV;
    xASL_adm_UnzipNifti(image.fname, 1); % unzip the NIfTI if needed
    spm_jobman('run', matlabbatch); % run the volumetrics
    
    % Edit the header names
    CSV = xASL_csvRead(SaveFileCSV);
    tableSPM12 = {'File' 'GM_volume_L' 'WM_volume_L' 'CSF_volume_L' 'Tissue_volume_L' 'Bone_volume_L' 'Air_volume_L'};
    tableSPM12(2,2:end) = CSV(2,size(CSV,2)-5:end);
    tableSPM12{2,1} = x.P.SubjectID;
    xASL_tsvWrite(tableSPM12, SaveFileTSV, 1);
    xASL_delete(SaveFileCSV);
else
    warning('No CAT12/SPM12 volumetric result file found, need to repeat segmentation');
end
if exist(SaveFileTSV, 'file')
    xASL_bids_csv2tsvReadWrite(SaveFileTSV);
end
%% 2) File management & store volumes for WMH segmentation (if exists)
if xASL_exist(x.P.Path_WMH_SEGM, 'file')
    xASL_adm_UnzipNifti(x.P.Path_WMH_SEGM); % make sure it is unzipped
    % This function prints the WMH volume (mL) and number of WMH in a TSV file
    bVerbose = 1; % verbal output
    BinThresh = 0.05; % default setting = 0.5. we did a cleanup procedure already, so can allow a very low threshold
    MinimalLesionVolume = 0.015; % default (mL), determines the minimal size of each detected lesion
    
    % Delete any previous volume files:
    xASL_adm_DeleteFileList(x.dir.SUBJECTDIR, '^LST_tlv_.*\.(csv|tsv)$', [], [0 Inf]);
    xASL_adm_DeleteFileList(x.D.TissueVolumeDir, '^LST_tlv_.*\.(csv|tsv)$', [], [0 Inf]);
    xASL_adm_DeleteFileList(x.D.TissueVolumeDir, ['WMH_LST_(LGA|LPA)_' x.P.SubjectID '.(csv|tsv)'], [], [0 Inf]);
    
    % Run the LST lesion volume calculator
	% -> Here, we change the current directory to the SUBJECT DIRECTORY as ps_LST_tlv outputs to the current (SUBJECT) directory
    CurrDir = pwd;
    cd(x.dir.SUBJECTDIR);
    [~, FileName] = ps_LST_tlv(x.P.Path_WMH_SEGM, ~bVerbose, BinThresh, MinimalLesionVolume);
    cd(CurrDir);
    
    % Get the current filename in the SUBJECT DIRECTORY
    FileName = fullfile(x.dir.SUBJECTDIR, FileName);
    % Define new filename
    pathOutput = fullfile(x.D.TissueVolumeDir, ['WMH_LST_' x.WMHsegmAlg '_' x.P.SubjectID '.tsv']);
    
    if ~exist(FileName, 'file')
         warning(['Missing:' FileName]);
    else
        tableLST = xASL_csvRead(FileName);
        tableLST{2,1} = x.P.SubjectID;
        xASL_tsvWrite(tableLST, pathOutput, true); % Write new output file
        xASL_delete(FileName);
    end
end
end