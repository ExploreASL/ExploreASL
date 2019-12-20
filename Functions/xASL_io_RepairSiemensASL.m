function xASL_io_RepairSiemensASL(InputPath)
%xASL_io_RepairSiemensASL Fix incorrect MOSAIC reordering of dcm2NIfTI slices/volumes for Siemens ASL
%
% FORMAT: xASL_io_RepairSiemensASL(AnalysisDir)
%
% INPUT:
%   AnalysisDir - path to folder containing the NIfTI/BIDS data (e.g. /data/RAD/share/EPAD/analysis) (REQUIRED)
%
% OUTPUT: n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function fixes incorrect reordering of dcm2NIfTI slices/volumes for ASL
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_io_RepairSiemensASL(AnalysisDir);
% __________________________________
% Copyright 2015-2019 ExploreASL



%% -------------------------------------------------------------------------------------
%% Admin
[Fdir, Ffile, Fext] = xASL_fileparts(InputPath);
BackupPath = fullfile(Fdir, [Ffile '_Backup' Fext]);

% Backup file
xASL_Copy(InputPath,BackupPath);

%% -------------------------------------------------------------------------------------
%% Run image repair
IM = xASL_io_Nifti2Im(InputPath);
NewIM = zeros(size(IM,1)*2,size(IM,2)*2,size(IM,3)/4,size(IM,4));

SzX = size(IM,1);
SzY = size(IM,2);

Sq1 = [1     SzX   SzY+1 SzY*2];
Sq2 = [SzX+1 SzX*2 SzY+1 SzY*2];
Sq3 = [1     SzX   1     SzY];
Sq4 = [SzX+1 SzX*2 1     SzY];

nSlices = 24;% size(IM,3)/4; % 36
% Assume that it runs in 24

for iSlice=1:size(IM,3)/nSlices
    TrackIndex = (iSlice-1)*24;
    Ind1 = TrackIndex+1;
    Ind2 = TrackIndex+2;
    Ind3 = TrackIndex+1+12;
    Ind4 = TrackIndex+2+12;

    Ind1b = TrackIndex+12-1;
    Ind2b = TrackIndex+12-0;
    Ind3b = TrackIndex+24-1;
    Ind4b = TrackIndex+24-0;

    IndA = (iSlice-1)*6+1;
    IndB = (iSlice-0)*6+0;

    % 1)
    NewIM(Sq1(1):Sq1(2),Sq1(3):Sq1(4),IndA:IndB ,:) = IM(:,:,[Ind1:2:Ind1b],:);

    % 2)
    NewIM(Sq2(1):Sq2(2),Sq2(3):Sq2(4),IndA:IndB,:) = IM(:,:,[Ind2:2:Ind2b],:);

    % 3)
    NewIM(Sq3(1):Sq3(2),Sq3(3):Sq3(4),IndA:IndB,:) = IM(:,:,[Ind3:2:Ind3b],:);

    % 4)
    NewIM(Sq4(1):Sq4(2),Sq4(3):Sq4(4),IndA:IndB,:) = IM(:,:,[Ind4:2:Ind4b],:);
end

%% -------------------------------------------------------------------------------------
%% Save
xASL_io_SaveNifti(InputPath,InputPath,NewIM,[],0);


end
