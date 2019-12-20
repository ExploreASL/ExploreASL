function xASL_io_RepairSiemensASL_MOSAIC(InputPath)
%xASL_io_RepairSiemensASL_MOSAIC Fix incorrect MOSAIC reordering of dcm2NIfTI slices/volumes for Siemens ASL
%
% FORMAT: xASL_io_RepairSiemensASL_MOSAIC(AnalysisDir)
%
% INPUT:
%   AnalysisDir - path to folder containing the NIfTI/BIDS data (e.g. /data/RAD/share/EPAD/analysis) (REQUIRED)
%
% OUTPUT: n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function fixes incorrect reordering of dcm2NIfTI MOSAIC slices/volumes for ASL
%              We need to know the original matrix size, which we try to
%              get with the AcquisitionMatrix field from the DICOM header
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_io_RepairSiemensASL_MOSAIC(AnalysisDir);
% __________________________________
% Copyright 2015-2019 ExploreASL



%% -------------------------------------------------------------------------------------
%% Admin
[Fdir, Ffile, Fext] = xASL_fileparts(InputPath);
% Get resolution from parms file, containing DICOM header AcquisitionMatrix
ParmsPath = xASL_adm_GetFileList(Fdir,'ASL4D.*_parms.mat', 'FPList', [0 Inf]); % do it this way, as sometimes the parms file had a different name, & AcquisitionMatrix should be same for all ASL files

if  isempty(ParmsPath)
    warning('Couldnt find parms file with AcquisitionMatrix for repairing Siemens MOSAIC:');
    fprintf('%s\n',InputPath);
    return;
else
    parms = load(ParmsPath{1}, '-mat');
    if ~isfield(parms.parms,'AcquisitionMatrix')
        warning('Couldnt find parms file with AcquisitionMatrix for repairing Siemens MOSAIC:');
        fprintf('%s\n',InputPath);
        return;
    else
        Resolution = 2*parms.parms.AcquisitionMatrix(1);
    end
end

% Backup file
BackupPath = fullfile(Fdir, [Ffile '_Backup' Fext]);
xASL_Copy(InputPath,BackupPath);


%% -------------------------------------------------------------------------------------
%% Run image repair
OldIM = xASL_io_Nifti2Im(InputPath);

nSlices = size(OldIM,1)^2/Resolution^2;
nRows = nSlices^0.5;
nColumns = nRows;

IsSquare2 = nRows/2 == round(nRows/2);
IsSquare = size(OldIM,1)==size(OldIM,2); % is each slice square?
Is1Slice = size(OldIM,3)==1; % are all slices contained in a single MOSAIC slice
IsEven = size(OldIM,1)/2==round(size(OldIM,1)/2); % are there an even number of voxels

if ~IsSquare2 || ~IsSquare || ~Is1Slice || ~IsEven % some basic checks
    warning('Couldnt repair this Siemens MOSAIC NifTI:');
    fprintf('%s\n',InputPath);
    return;
end

OldIM = xASL_im_rotate(OldIM, 90);

for iVolume=1:size(OldIM,4)
    iSlice=1;
    for iRow=1:nRows
        for iColumn=1:nColumns
            xStart = (iColumn-1)*Resolution+1;
            xEnd = iColumn*Resolution;
            yStart = (iRow-1)*Resolution+1;
            yEnd = iRow*Resolution;
            NewIM(:,:,iSlice,iVolume) = OldIM(yStart:yEnd,xStart:xEnd,1,iVolume);
            iSlice = iSlice+1;
        end
    end
end


%% -------------------------------------------------------------------------------------
%% Save
xASL_io_SaveNifti(InputPath, InputPath, NewIM, [], false);


end
