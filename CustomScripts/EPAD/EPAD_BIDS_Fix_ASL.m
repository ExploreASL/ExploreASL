function EPAD_BIDS_Fix_ASL(AnalysisDir)
%EPAD_BIDS_Fix_ASL Fix incorrect reordering of dcm2NIfTI slices/volumes for ASL
%
% FORMAT: EPAD_BIDS_Fix_ASL(AnalysisDir)
% 
% INPUT:
%   AnalysisDir - path to folder containing the NIfTI/BIDS data (e.g. /data/RAD/share/EPAD/analysis) (REQUIRED)
%
% OUTPUT: n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function fixes incorrect reordering of dcm2NIfTI slices/volumes for ASL
%              Currently for site 110, 012 and 031
%              110 provided MOSAIC DICOMs which convert poorly with
%              dcm2niiX, which we correct here
%              012 & 031 have a first volume that should be removed
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: EPAD_BIDS_Fix_ASL(AnalysisDir);
% __________________________________
% Copyright 2015-2019 ExploreASL

%% ASL site 10 has wrong slice ordering, fix this
fprintf('%s','Fixing ASL slice ordering for site 010:   ');
FileList = xASL_adm_GetFsList(AnalysisDir, '^(010|110)-\d*$', 1, [], [], [0 Inf]);

for iFile=1:length(FileList)
    xASL_TrackProgress(iFile, length(FileList));
    ASLfile     = fullfile(AnalysisDir, FileList{iFile}, 'ASL_1','ASL4D.nii');
    ASLBackup   = fullfile(AnalysisDir, FileList{iFile}, 'ASL_1','ASL4D_Backup.nii');
    % Check if Backup file already exists, that means this file has already been done
    if  xASL_exist(ASLfile, 'file') && ~ xASL_exist(ASLBackup, 'file')
        xASL_io_RepairSiemensASL_MOSAIC(ASLfile);
    end
end
fprintf('\n');

%% First image of Siemens ASL site 12/31 should be removed
%  Which we do here by splitting it into an M0 and deleting this M0 NIfTI
%  afterwards
fprintf('%s','Fixing ASL slice ordering for site 012/031:   ');
FileList   = xASL_adm_GetFsList(AnalysisDir, '^(012|031)-\d{4}$', 1, [], [], [0 Inf]);

for iFile=1:length(FileList)
    xASL_TrackProgress(iFile, length(FileList));
    
    PathASL         = fullfile(AnalysisDir, FileList{iFile}, 'ASL_1','ASL4D.nii');
    PathASLBackup   = fullfile(AnalysisDir, FileList{iFile}, 'ASL_1','ASL4D_Backup.nii');
    PathM0          = fullfile(AnalysisDir, FileList{iFile}, 'ASL_1','M0.nii');
    PathM0parms     = fullfile(AnalysisDir, FileList{iFile}, 'ASL_1','M0_parms.mat');    

    if xASL_exist(PathASL, 'file') && ~ xASL_exist(PathASLBackup, 'file')
        xASL_Copy(PathASL, PathASLBackup);
        xASL_delete(PathM0);
        xASL_delete(PathM0parms);

        xASL_io_SplitASL_M0(PathASL, 1);
    end
end
fprintf('\n');


end

